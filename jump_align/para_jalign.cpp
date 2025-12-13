#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <climits>
#include <parasail.h>
#include "logging.h"

using namespace std;

// some hard limits
#define QNAME_MAX_LEN 1000
#define SEQ_MAX_LEN	1000000
#define REF_MAX_LEN 1000000

// calculated limits
#define LINE_MAX_LEN (QNAME_MAX_LEN + 1 + SEQ_MAX_LEN + 1 + REF_MAX_LEN + 1 + REF_MAX_LEN + 1 + 1)

// parametized limits
#define MIN_ACCEPTABLE_JREADLEN 30

// line buffers used for reading input
// double buffering used to optimize the (common?) case where the reference repeats from the previous line
static char linebuf[2][LINE_MAX_LEN];
static int linebuf_index = 0;

// handy structure to keep reference pointers and lengths
typedef struct {
	char* ptr[2];
	int len[2];
} refs_t;

int s_match;
int s_mismatch;
int s_open;
int s_gap; 
int s_jump_penalty;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <parasail.h>
#include <parasail/matrix_lookup.h>

typedef struct {
    int score;
    int readlen1;
    int readlen2;
    int score1;
    int score2;
} result_t;

// Reverse a string in-place
static void reverse_string(char *s) {
    int len = strlen(s);
    for (int i = 0; i < len / 2; i++) {
        char tmp = s[i];
        s[i] = s[len - 1 - i];
        s[len - 1 - i] = tmp;
    }
}

// Compute best score per query position (F vector)
static int *compute_F(const char *query, int m, const char *ref, int n, const parasail_matrix_t *matrix) {
    parasail_result_t *result = parasail_sw_table(query, m, ref, n, 10, 1, matrix);
    const int *table = parasail_result_get_score_table(result);
    int *F = (int *)calloc(m, sizeof(int));

    for (int i = 0; i < m; i++) {
        int max_row = 0;
        for (int j = 0; j < n; j++) {
            int val = table[i*n + j];
            if (val > max_row) max_row = val;
        }
        F[i] = max_row;
    }

    parasail_result_free(result);
    return F;
}

result_t jalign(const char *query, const char *ref1, const char *ref2) {

    result_t result;

    int nquery = strlen(query);
    int nref1 = strlen(ref1);
    int nref2 = strlen(ref2);

#ifdef DEBUG_PRINTS
    printf("\n--- query (len %d) ---\n", nquery);
    printf("%s\n", query);
    printf("\n--- ref1 (len %d) ---\n", nref1);
    printf("%s\n", ref1);
    printf("\n--- ref2 (len %d) ---\n", nref2);
    printf("%s\n", ref2);
#endif    

    parasail_matrix_t *matrix = parasail_matrix_create("ACGT", s_match, s_mismatch);

    // --- First alignment ---
    parasail_result_t *result1 = parasail_sw_table_diag_16(query, nquery, ref1, nref1, -s_open, -s_gap, matrix);
#ifdef DEBUG_PRINTS
    parasail_result_t *result1_trace = parasail_sw_trace_diag_16(query, m, ref1, n1, -s_open, -s_gap, matrix);
#endif
    const int *table1 = parasail_result_get_score_table(result1);
    int *forward = (int *)calloc(nquery, sizeof(int));
    const int *pcell = table1;
    int cell_val;
    for (int row = 0; row < nquery; row++) {
        int max_val = INT_MIN;
        for (int j = 0; j < nref1; j++) {
            if ( (cell_val = *pcell++) > max_val )
                max_val = cell_val;
        }
        forward[row] = max_val;
    }

    // --- Second alignment (reverse phase) ---
    char *rev_query = (char *)strdup(query);
    reverse_string(rev_query);
    char *rev_ref2 = (char *)strdup(ref2);
    reverse_string(rev_ref2);

    parasail_result_t *result2 = parasail_sw_table_diag_16(rev_query, nquery, rev_ref2, nref2, -s_open, -s_gap, matrix);
#ifdef DEBUG_PRINTS
    parasail_result_t *result2_trace = parasail_sw_trace_diag_16(rev_query, m, rev_ref2, n2, -s_open, -s_gap, matrix);
#endif
    const int *table2 = parasail_result_get_score_table(result2);
    int *backward = (int *)calloc(nquery, sizeof(int));
    pcell = table2;
    for (int i = 0; i < nquery; i++) {
        int max_val = INT_MIN;
        for (int j = 0; j < nref2; j++) {
            if ( (cell_val = *pcell++) > max_val )
                max_val = cell_val;
        }
        backward[nquery - i - 1] = max_val; // reverse back to original order
    }

    // --- Combine ---
    int best_score = INT_MIN;
    int best_q = -1;
    for (int q = MIN_ACCEPTABLE_JREADLEN ; q < nquery - 1 - MIN_ACCEPTABLE_JREADLEN ; q++) {
        int s = forward[q] + backward[q + 1] - s_jump_penalty;
        if (s > best_score) {
            best_score = s;
            best_q = q;
        }
    }

    result.score = best_score;
    result.score1 = result1->score;
    result.score2 = result2->score;
    result.readlen1 = best_q;
    result.readlen2 = nquery - best_q;

    #ifdef DEBUG_PRINTS
    printf("Best jump at query position %d, combined score = %d, score1 %d sscore2 %d better %d\n", 
        best_q + 1, best_score, result1->score, result2->score, better);
        
    // --- Tracebacks ---
    parasail_traceback_t *tb1 = parasail_result_get_traceback(
        result1_trace, query, nquery, ref1, nref1, matrix, '|', '+', ' ');

    parasail_traceback_t *tb2 = parasail_result_get_traceback(
        result2_trace, rev_query, nquerym, rev_ref2, nref2, matrix, '|', '+', ' ');

    // Reverse traceback strings for the second alignment
    reverse_string(tb2->query);
    reverse_string(tb2->comp);
    reverse_string(tb2->ref);

    printf("\n--- Alignment 1 (pre-jump, len %ld %ld %ld) ---\n",
                strlen(tb1->query), strlen(tb1->comp), strlen(tb1->ref));
    printf("%s\n%s\n%s\n", tb1->query, tb1->comp, tb1->ref);

    printf("\n--- Alignment 2 (post-jump, len %ld %ld %ld) ---\n",
                strlen(tb2->query), strlen(tb2->comp), strlen(tb2->ref));
    printf("%s\n%s\n%s\n", tb2->query, tb2->comp, tb2->ref);
#endif

    // --- Cleanup ---
    free(forward);
    free(backward);
    free(rev_query);
    free(rev_ref2);
    parasail_result_free(result1);
    parasail_result_free(result2);
#ifdef DEBUG_PRINTS
    parasail_traceback_free(tb1);
    parasail_traceback_free(tb2);
    parasail_result_free(result1_trace);
    parasail_result_free(result2_trace);
#endif
    parasail_matrix_free(matrix);

    return result;
}

// main entry point. read input csv (SEQ REF1 REF2) and generate alignment output CSV
int main(int argc, char* argv[]) {
	
	// check args
	if ( argc < 7 ) {
		fprintf(stderr, "wrong score arguments:\n");
		fprintf(stderr, "usage: %s <match> <mismatch> <open> <extend> <do_del> <jump> [input_file]\n", argv[0]);
		exit(-1);
	}

    // extract scores
    s_match = atoi(argv[1]);
    s_mismatch = atoi(argv[2]);
    s_open = atoi(argv[3]);
    s_gap = atoi(argv[4]);
    s_jump_penalty = atoi(argv[6]);

	// open input file
	FILE* infile;
	if ( argc >= 8 ) {
		if ( !(infile = fopen(argv[7], "r")) ) {
			fprintf(stderr, "failed to open input file: %s\n", argv[7]);
			exit(-1);
		}
	} else {
		infile = stdin;
	}


	// loop on reading lines and process
	int lineno = 0;
	refs_t last_refs;
	printf("readName\t");
	printf("score\tjreadlen1\tjreadlen2\t");
	printf("dscore\tdjreadlen1\tdjreadlen2\t");
	printf("score1\tscore2\n");
	while ( fgets(linebuf[linebuf_index], sizeof(linebuf[0]), infile) ) {
		lineno++;
		char* line = linebuf[linebuf_index];

		// ignore empty or comment lines
		if ( !line[0] || line[0] == '\n' || line[0] == '#' ) {
			continue;
		}

		// ignore lines not starting with a alphanum
		if ( !isalnum(line[0]) ) {
			continue;
		}

		// parse the line, get qname and sequence
		char* qname = strtok(line, "\t");
		if ( !qname ) {
			LOG(ERR) << "failed to read qname on lineno: " << lineno;
		}
		char* seq = strtok(nullptr, "\t");
		if ( !seq ) {
			LOG(ERR) << "failed to read sequence on lineno: " << lineno;
		}
		int seq_len = strlen(seq);
		if ( seq_len > SEQ_MAX_LEN ) {
			LOG(WARN) << "sequence is longer than defined max length: " << seq_len << " > " << SEQ_MAX_LEN;
		}

		// get references
		bool last_refs_used = false;
		refs_t refs;
		for ( int i = 0 ; i < 2 && !last_refs_used ; i++ ) {
			if ( !(refs.ptr[i] = strtok(nullptr, "\t\n")) ) {
				LOG(ERR) << "failed to reference no " << i << " on lineno: " << lineno;
			}
			if ( refs.ptr[i][0] == '=' ) {
				refs = last_refs;
				last_refs_used = true;
			} else {
				if ( (refs.len[i] = strlen(refs.ptr[i])) > REF_MAX_LEN ) {
					LOG(WARN) << "reference no " << i << " is longer than defined max length: " << refs.len[i] << " > " << REF_MAX_LEN;
				}
			}
		}

        // jump align
        ///
        result_t result = jalign(seq, refs.ptr[0], refs.ptr[1]);
        result_t dresult = jalign(seq, refs.ptr[1], refs.ptr[0]);

		// process jump align results
		printf("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", 
            qname, 
            result.score, result.readlen1, result.readlen2,
            dresult.score, dresult.readlen1, dresult.readlen2,
            result.score1, result.score2);

		// save last references, switch buffers if did not use last
		last_refs = refs;		
		if ( !last_refs_used ) {
			linebuf_index = 1 - linebuf_index;
		}
	}

	return 0;

}

#ifdef TEST_CODE
typedef struct {
    int start1, end1;
    int start2, end2;
    int score;
} alignment_t;

typedef struct {
    alignment_t first;
    alignment_t second;
    int total_score;
} jump_alignment_t;

static alignment_t find_best_local(
    const char* s1, const char* s2,
    const parasail_matrix_t* matrix,
    int open, int extend,
    int* mask1, int* mask2)
{
    int n = strlen(s1);
    int m = strlen(s2);

    // Build masked copies
    char* masked1 = (char*)malloc(n + 1);
    char* masked2 = (char*)malloc(m + 1);
    for (int i = 0; i < n; i++) masked1[i] = mask1 && mask1[i] ? 'X' : s1[i];
    for (int j = 0; j < m; j++) masked2[j] = mask2 && mask2[j] ? 'X' : s2[j];
    masked1[n] = masked2[m] = '\0';

    parasail_result_t* result = parasail_sw(masked1, n, masked2, m, open, extend, matrix);

    alignment_t aln = {0};
    aln.score = parasail_result_get_score(result);
    aln.end1  = parasail_result_get_end_query(result);
    aln.end2  = parasail_result_get_end_ref(result);

    // Estimate starts via traceback
    parasail_result_t* trace = parasail_sw_trace(masked1, n, masked2, m, open, extend, matrix);
    parasail_cigar_t* cigar = parasail_result_get_cigar(trace, masked1, n, masked2, m, matrix);
    int len = cigar->len;
    int pos1 = aln.end1;
    int pos2 = aln.end2;

    for (int k = len - 1; k >= 0; k--) {
        int op_len = cigar->seq[k] >> 4;
        char op = cigar->seq[k] & 15;
        if (op == 0 || op == 7 || op == 8) { pos1 -= op_len; pos2 -= op_len; }
        else if (op == 1) pos2 -= op_len;
        else if (op == 2) pos1 -= op_len;
        if (pos1 <= 0 || pos2 <= 0) break;
    }
    aln.start1 = pos1;
    aln.start2 = pos2;

    // Mask region
    if (mask1)
        for (int i = aln.start1; i <= aln.end1 && i < n; i++)
            mask1[i] = 1;
    if (mask2)
        for (int j = aln.start2; j <= aln.end2 && j < m; j++)
            mask2[j] = 1;

    parasail_cigar_free(cigar);
    parasail_result_free(trace);
    parasail_result_free(result);
    free(masked1);
    free(masked2);
    return aln;
}

jump_alignment_t jump_local_align_parasail(
    const char* seq1,
    const char* seq2,
    const parasail_matrix_t* matrix,
    int open, int extend)
{
    int n = strlen(seq1), m = strlen(seq2);
    int* mask1 = (int*)calloc(n, sizeof(int));
    int* mask2 = (int*)calloc(m, sizeof(int));

    alignment_t a1 = find_best_local(seq1, seq2, matrix, open, extend, mask1, mask2);
    alignment_t a2 = find_best_local(seq1, seq2, matrix, open, extend, mask1, mask2);

    jump_alignment_t out = { a1, a2, a1.score + a2.score };
    free(mask1);
    free(mask2);
    return out;
}

int main(void)
{
    const char* seq1 = "ACGTACGTGAC";
    const char* seq2 = "TACGTGTTACGTGAC";

    parasail_matrix_t* matrix = parasail_matrix_create("ACGT", 2, -1);
    jump_alignment_t j = jump_local_align_parasail(seq1, seq2, matrix, 10, 1);

    printf("Segment 1: seq1[%d..%d] ↔ seq2[%d..%d], score=%d\n",
           j.first.start1, j.first.end1, j.first.start2, j.first.end2, j.first.score);
    printf("Segment 2: seq1[%d..%d] ↔ seq2[%d..%d], score=%d\n",
           j.second.start1, j.second.end1, j.second.start2, j.second.end2, j.second.score);
    printf("Total score: %d\n", j.total_score);

    parasail_matrix_free(matrix);
    return 0;
}

#endif
