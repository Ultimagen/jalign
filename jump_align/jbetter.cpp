#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>
#include <htslib/hts_log.h>
#include <parasail.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <climits>
#include <cstdlib>
#include <cctype>

#define REF_PAD 500
#define MIN_MISMATCHES 20

typedef enum {
    DEL,
    DUP
} jmode_t;

jmode_t mode = DUP;

int s_match = 10;
int s_mismatch = -25;
int s_open = -10;
int s_gap = -5;
int s_jump_penalty = 0;

typedef struct {
    int jscore;
    int score1;
    int score2;
    int jumpAt;
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

result_t jalign_better(const char *query, const char *ref1, const char *ref2, char* ref2_rev) {

    result_t result;

    int m = strlen(query);
    int n1 = strlen(ref1);
    int n2 = strlen(ref2);

#ifdef DEBUG_PRINTS
    printf("\n--- query (len %d) ---\n", m);
    printf("%s\n", query);
    printf("\n--- ref1 (len %d) ---\n", n1);
    printf("%s\n", ref1);
    printf("\n--- ref2 (len %d) ---\n", n2);
    printf("%s\n", ref2);
#endif    

    parasail_matrix_t *matrix = parasail_matrix_create("ACGT", s_match, s_mismatch);

    // --- First alignment ---
    parasail_result_t *result1 = parasail_sw_table_diag_16(query, m, ref1, n1, -s_open, -s_gap, matrix);
    assert(result1);
#ifdef DEBUG_PRINTS
    parasail_result_t *result1_trace = parasail_sw_trace_diag_16(query, m, ref1, n1, -s_open, -s_gap, matrix);
#endif
    const int *table1 = parasail_result_get_score_table(result1);
    int *F = (int *)calloc(m, sizeof(int));
    for (int i = 0; i < m; i++) {
        int max_row = 0;
        for (int j = 0; j < n1; j++) {
            int val = table1[i*n1 + j];
            if (val > max_row) max_row = val;
        }
        F[i] = max_row;
    }

    // --- Second alignment (reverse phase) ---
    char *rev_query = (char *)strdup(query);
    reverse_string(rev_query);
    char *rev_ref2;
    if ( ref2_rev ) {
        rev_ref2 = ref2_rev;
    }
    else {
        rev_ref2 = (char *)strdup(ref2);
        reverse_string(rev_ref2);
    }

    parasail_result_t *result2 = parasail_sw_table_diag_16(rev_query, m, rev_ref2, n2, -s_open, -s_gap, matrix);
    assert(result2);
#ifdef DEBUG_PRINTS
    parasail_result_t *result2_trace = parasail_sw_trace_diag_16(rev_query, m, rev_ref2, n2, -s_open, -s_gap, matrix);
#endif
    const int *table2 = parasail_result_get_score_table(result2);
    int *B = (int *)calloc(m, sizeof(int));
    for (int i = 0; i < m; i++) {
        int max_row = 0;
        for (int j = 0; j < n2; j++) {
            int val = table2[i*n2 + j];
            if (val > max_row) max_row = val;
        }
        B[m - i - 1] = max_row; // reverse back to original order
    }

    // --- Combine ---
    int best_score = 0;
    int best_q = -1;
    for (int q = 0; q < m - 1; q++) {
        int s = F[q] + B[q + 1] - s_jump_penalty;
        if (s > best_score) {
            best_score = s;
            best_q = q;
        }
    }

    result.jscore = best_score;
    result.score1 = result1->score;
    result.score2 = result2->score;
    result.jumpAt = best_q + 1;

    #ifdef DEBUG_PRINTS
    printf("Best jump at query position %d, combined score = %d, score1 %d sscore2 %d better %d\n", 
        best_q + 1, best_score, result1->score, result2->score, better);
        
    // --- Tracebacks ---
    parasail_traceback_t *tb1 = parasail_result_get_traceback(
        result1_trace, query, m, ref1, n1, matrix, '|', '+', ' ');

    parasail_traceback_t *tb2 = parasail_result_get_traceback(
        result2_trace, rev_query, m, rev_ref2, n2, matrix, '|', '+', ' ');

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
    free(F);
    free(B);
    free(rev_query);
    if ( !ref2_rev ) {
        free(rev_ref2);
    }
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

// Forward declaration of the read-level function
bool jbetter_read(const std::string& read_seq,
                  const std::string& ref_start_seq,
                  const std::string& ref_end_seq) {
    const std::string& ref1 = (mode == DUP) ? ref_end_seq : ref_start_seq;
    const std::string& ref2 = (mode == DUP) ? ref_start_seq : ref_end_seq;
    result_t r = jalign_better(read_seq.c_str(), ref1.c_str(), ref2.c_str(), nullptr);
    return r.jscore > std::max(r.score1, r.score2);
}

int jbetter(const std::vector<bam1_t*>& reads,
            const std::string& ref_start_seq,
            const std::string& ref_end_seq)
{
    int score = 0;

    for (const bam1_t* b : reads) {
        const uint8_t* seq_data = bam_get_seq(b);
        int len = b->core.l_qseq;

        std::string read_seq;
        read_seq.reserve(len);

        for (int i = 0; i < len; ++i) {
            char base = seq_nt16_str[bam_seqi(seq_data, i)];
            read_seq.push_back(base);
        }

        assert(ref_start_seq.size());
        assert(ref_end_seq.size());
        if (jbetter_read(read_seq, ref_start_seq, ref_end_seq))
            ++score;
    }

    return score;
}

int gs_init() {
    FILE *fp;
    char token[1024];

    // Run the command and open a pipe to read its output
    fp = popen("gcloud auth print-access-token", "r");
    if (fp == NULL) {
        perror("popen failed");
        return -1;
    }

    // Read the first line of output (the token)
    if (fgets(token, sizeof(token), fp) == NULL) {
        perror("fgets failed");
        pclose(fp);
        return -1;
    }

    // Close the pipe
    if (pclose(fp) == -1) {
        perror("pclose failed");
        return -1;
    }

    // Remove trailing newline if present
    size_t len = strlen(token);
    if (len > 0 && token[len - 1] == '\n') {
        token[len - 1] = '\0';
    }

    // Set the environment variable
    if (setenv("GCS_OAUTH_TOKEN", token, 1) != 0) {
        perror("setenv failed");
        return -1;
    }

    return 0;
}

// BED region structure
struct BedRegion {
    std::string chrom;
    int start;
    int end;
    std::string rest;
};

// Parse one BED line
bool parse_bed_line(const std::string& line, BedRegion& region) {
    if (line.empty() || line[0] == '#') return false;
    std::istringstream ss(line);
    ss >> region.chrom >> region.start >> region.end;
    std::getline(ss, region.rest);
    if (!region.rest.empty() && region.rest[0] == '\t')
        region.rest.erase(0, 1);
    return true;
}

bool read_accepted(bam1_t* b, faidx_t* fai, const bam_hdr_t* hdr, int N)
{
    const uint32_t *cigar = bam_get_cigar(b);
    const uint8_t *seq = bam_get_seq(b);
    int n_cigar = b->core.n_cigar;

    int softclip_bases = 0;
    for (int i = 0; i < n_cigar; ++i) {
        int op  = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CSOFT_CLIP)
            softclip_bases += len;
    }

    // Early exit if already exceeds threshold
    if (softclip_bases >= N)
        return true;

    // If unmapped, cannot compare to reference
    if (b->core.tid < 0)
        return false;

    const char* ref_name = hdr->target_name[b->core.tid];
    int ref_start = b->core.pos;  // 0-based
    int ref_end = bam_endpos(b) - 1;

    int ref_len = 0;
    char* ref_seq = faidx_fetch_seq(fai, ref_name, ref_start, ref_end, &ref_len);
    if (!ref_seq)
        return false;  // failed to get reference

    // Compare aligned positions
    int read_pos = 0, ref_pos = 0;
    int mismatch_bases = 0;

    for (int i = 0; i < n_cigar; ++i) {
        int op  = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        switch (op) {
        case BAM_CMATCH:
        case BAM_CEQUAL:
        case BAM_CDIFF:
            for (int j = 0; j < len; ++j) {
                if (ref_pos + j >= ref_len) break;

                char read_base = seq_nt16_str[bam_seqi(seq, read_pos + j)];
                char ref_base  = std::toupper(ref_seq[ref_pos + j]);
                if (read_base != ref_base && read_base != 'N' && ref_base != 'N')
                    mismatch_bases++;
            }
            read_pos += len;
            ref_pos  += len;
            break;

        case BAM_CINS:
            // insertion relative to ref → all inserted bases count as mismatches
            mismatch_bases += len;
            read_pos += len;
            break;

        case BAM_CDEL:
        case BAM_CREF_SKIP:
            ref_pos += len;
            break;

        case BAM_CSOFT_CLIP:
            // already counted above
            read_pos += len;
            break;

        case BAM_CHARD_CLIP:
        case BAM_CPAD:
        default:
            break;
        }
    }

    free(ref_seq);

    int total = softclip_bases + mismatch_bases;
    return (total >= N);
}

// Collect reads overlapping a specific coordinate
std::vector<bam1_t*> fetch_reads_at(hts_idx_t* idx, samFile* fp, bam_hdr_t* hdr,
                                    int tid, int pos, faidx_t* fai)
{
    std::vector<bam1_t*> reads;

    hts_itr_t* iter = sam_itr_queryi(idx, tid, pos, pos + 1);
    if (!iter) return reads;

    bam1_t* b = bam_init1();
    while (sam_itr_next(fp, iter, b) >= 0) {
        int read_start = b->core.pos;
        int read_end = bam_endpos(b);
        if ( !read_accepted(b, fai, hdr, MIN_MISMATCHES) ) {
            continue;
        }
        if (pos >= read_start && pos < read_end) {
            bam1_t* copy = bam_dup1(b);
            reads.push_back(copy);
        }
    }

    bam_destroy1(b);
    hts_itr_destroy(iter);
    return reads;
}

// Compute reference span (min start / max end) from reads
std::pair<int,int> compute_span(const std::vector<bam1_t*>& reads, int initial_center) {
    int min_start = initial_center - REF_PAD, max_end = initial_center + REF_PAD;
    for (auto b : reads) {
        min_start = std::min(min_start, (int)b->core.pos);
        max_end = std::max(max_end, (int)bam_endpos(b));
    }
    return {min_start, max_end};
}

// Free bam1_t* vector
void destroy_reads(std::vector<bam1_t*>& reads) {
    for (auto b : reads) bam_destroy1(b);
    reads.clear();
}

// Deduplicate reads by QNAME
std::vector<bam1_t*> deduplicate_reads(const std::vector<bam1_t*>& reads) {
    std::unordered_set<std::string> seen;
    std::vector<bam1_t*> unique;
    for (auto b : reads) {
        std::string name(bam_get_qname(b));
        if (seen.insert(name).second) {
            unique.push_back(b);
        } else {
            // duplicate: free memory immediately
            bam_destroy1(b);
        }
    }
    return unique;
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " regions.bed reads.cram ref.fasta\n";
        return 1;
    }

    std::string bed_path = argv[1];
    std::string cram_path = argv[2];
    std::string ref_path = argv[3];

    // Open CRAM
    gs_init();
    samFile* fp = sam_open(cram_path.c_str(), "r");
    if (!fp) {
        std::cerr << "Error: failed to open CRAM " << cram_path << "\n";
        return 1;
    }

    bam_hdr_t* hdr = sam_hdr_read(fp);
    if (!hdr) {
        std::cerr << "Error: failed to read CRAM header.\n";
        sam_close(fp);
        return 1;
    }

    // Attach reference to CRAM
    //hts_set_log_level(HTS_LOG_DEBUG);
    if (hts_set_fai_filename(fp, ref_path.c_str()) != 0) {
        std::cerr << "Error: failed to attach reference FASTA.\n";
        bam_hdr_destroy(hdr);
        sam_close(fp);
        return 1;
    }

    // Load CRAM index
    hts_idx_t* idx = sam_index_load(fp, cram_path.c_str());
    if (!idx) {
        std::cerr << "Error: failed to load CRAM index (.crai).\n";
        bam_hdr_destroy(hdr);
        sam_close(fp);
        return 1;
    }

    // Load FASTA index
    faidx_t* fai = fai_load(ref_path.c_str());
    if (!fai) {
        std::cerr << "Error: failed to load FASTA index (.fai).\n";
        hts_idx_destroy(idx);
        bam_hdr_destroy(hdr);
        sam_close(fp);
        return 1;
    }

    std::ifstream bed(bed_path);
    if (!bed) {
        std::cerr << "Error: failed to open BED file " << bed_path << "\n";
        fai_destroy(fai);
        hts_idx_destroy(idx);
        bam_hdr_destroy(hdr);
        sam_close(fp);
        return 1;
    }

    std::string line;
    while (std::getline(bed, line)) {
        BedRegion region;
        if (!parse_bed_line(line, region)) continue;

        int tid = bam_name2id(hdr, region.chrom.c_str());
        if (tid < 0) {
            std::cerr << "Warning: chromosome " << region.chrom << " not found in CRAM.\n";
            continue;
        }

        // Fetch reads crossing region start and end
        auto start_reads = fetch_reads_at(idx, fp, hdr, tid, region.start, fai);
        auto end_reads   = fetch_reads_at(idx, fp, hdr, tid, region.end - 1, fai);

        // Combine and deduplicate
        auto start_span = compute_span(start_reads, region.start);
        auto end_span   = compute_span(end_reads, region.end);
        std::vector<bam1_t*> all_reads = start_reads;
        all_reads.insert(all_reads.end(), end_reads.begin(), end_reads.end());
        std::vector<bam1_t*> unique_reads = deduplicate_reads(all_reads);


        // Fetch reference sequences
        int seqlen = 0;
        std::string start_seq, end_seq;
        if (start_span.first >= 0 && start_span.second > start_span.first) {
            char* seq = faidx_fetch_seq(fai, region.chrom.c_str(),
                                        start_span.first,
                                        start_span.second - 1,
                                        &seqlen);
            if (seq) {
                assert(strlen(seq));
                assert(seqlen);
                start_seq.assign(seq, seqlen);
                free(seq);
                assert(start_seq.size());
            }
        }
        if (end_span.first >= 0 && end_span.second > end_span.first) {
            char* seq = faidx_fetch_seq(fai, region.chrom.c_str(),
                                        end_span.first,
                                        end_span.second - 1,
                                        &seqlen);
            if (seq) {
                assert(strlen(seq));
                assert(seqlen);
                end_seq.assign(seq, seqlen);
                free(seq);
                assert(end_seq.size());
            }
        }

        // Call jbetter
        assert(start_seq.size());
        assert(end_seq.size());
        int score = jbetter(unique_reads, start_seq, end_seq);

        // Output BED line with result
        std::cout << region.chrom << '\t' << region.start << '\t' << region.end;
        if (!region.rest.empty()) std::cout << '\t' << region.rest;
        std::cout << '\t' << score << '\n';

        // Cleanup
        /*
        destroy_reads(start_reads);
        destroy_reads(end_reads);
        */
        destroy_reads(unique_reads);
    }

    fai_destroy(fai);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(fp);
    return 0;
}
