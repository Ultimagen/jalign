#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <parasail.h>

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
