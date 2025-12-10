#include <stdio.h>
#include <string.h>
#include <string>
#include "GlobalAligner.hpp"
#include "GlobalJumpAligner.hpp"
#include "logging.h"

// some hard limits
#define QNAME_MAX_LEN 1000
#define SEQ_MAX_LEN	1000000
#define REF_MAX_LEN 1000000

// calculated limits
#define LINE_MAX_LEN (QNAME_MAX_LEN + 1 + SEQ_MAX_LEN + 1 + REF_MAX_LEN + 1 + REF_MAX_LEN + 1 + 1)

// line buffers used for reading input
// double buffering used to optimize the (common?) case where the reference repeats from the previous line
static char linebuf[2][LINE_MAX_LEN];
static int linebuf_index = 0;

// handy structure to keep reference pointers and lengths
struct refs_t {
	char* ptr[2];
	int len[2];
};

using namespace std;

// convert an alignment path to a cigar string
string to_string(ALIGNPATH::path_t& path) {
	string s;
	for ( auto&& segment : path ) {
		s.append(to_string(segment.length));
		s += segment_type_to_cigar_code(segment.type);
	}
	return s;
}

// main entry point. read input csv (SEQ REF1 REF2) and generate alignment output CSV
int main(int argc, char* argv[]) {
	
	// check args
	if ( argc < 7 ) {
		fprintf(stderr, "wrong score arguments:\n");
		fprintf(stderr, "usage: %s <match> <mismatch> <open> <extend> <offedge> <jump> [input_file]\n", argv[0]);
		exit(-1);
	}

	// configure aligners
	AlignmentScores<int> scores(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), false);
	GlobalJumpAligner<int> jump_aligner(scores, atoi(argv[6]));
	GlobalAligner<int> aligner1(scores);
	GlobalAligner<int> aligner2(scores);

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
	struct refs_t last_refs;
	printf("readName\t");

	printf("score\tjumpInsertSize\tjumpRange\t");
	printf("jbegin1\tjapath1\tjreadlen1\tjreflen1\t");
	printf("jbegin2\tjapath2\tjreadlen2\tjreflen2\t");

	printf("dscore\tdjumpInsertSize\tdjumpRange\t");
	printf("djbegin1\tdjapath1\tdjreadlen1\tdjreflen1\t");
	printf("djbegin2\tdjapath2\tdjreadlen2\tdjreflen2\t");

	printf("score1\tbegin1\tapath1\treadlen1\treflen1\t");
	printf("score2\tbegin2\tapath2\treadlen2\treflen2\n");
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
		struct refs_t refs;
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

		// setup results for the different aligners
		JumpAlignmentResult<int> result;
		JumpAlignmentResult<int> dresult;
		AlignmentResult<int> result1;
		AlignmentResult<int> result2;


		jumpAlign<GlobalJumpAligner<int>,const char*,int>(
			jump_aligner, seq, seq + seq_len, refs.ptr[0], refs.ptr[0] + refs.len[0], refs.ptr[1], refs.ptr[1] + refs.len[1], result);
		jumpAlign<GlobalJumpAligner<int>,const char*,int>(
			jump_aligner, seq, seq + seq_len, refs.ptr[1], refs.ptr[1] + refs.len[1], refs.ptr[0], refs.ptr[0] + refs.len[0], dresult);
		aligner1.align(seq, seq + seq_len, refs.ptr[0], refs.ptr[0] + refs.len[0], result1);
		aligner2.align(seq, seq + seq_len, refs.ptr[1], refs.ptr[1] + refs.len[1], result2);

		// process jump align results
		string apath1 = to_string(result.align1.apath);
		string apath2 = to_string(result.align2.apath);
		printf("%s\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%s\t%d\t%d", 
					qname,
					result.score, result.jumpInsertSize, result.jumpRange,
					result.align1.beginPos, apath1.c_str(), apath_read_length(result.align1.apath), apath_ref_length(result.align1.apath),
					result.align2.beginPos, apath2.c_str(), apath_read_length(result.align2.apath), apath_ref_length(result.align2.apath)
					);

		// process second jump align results
		string dapath1 = to_string(dresult.align1.apath);
		string dapath2 = to_string(dresult.align2.apath);
		printf("\t%d\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%s\t%d\t%d", 
					dresult.score, dresult.jumpInsertSize, dresult.jumpRange,
					dresult.align1.beginPos, dapath1.c_str(), apath_read_length(dresult.align1.apath), apath_ref_length(dresult.align1.apath),
					dresult.align2.beginPos, dapath2.c_str(), apath_read_length(dresult.align2.apath), apath_ref_length(dresult.align2.apath)
					);


		// simple align results
		for ( int i = 0 ; i < 2 ; i++ ) {
			AlignmentResult<int>& result = !i ? result1 : result2;
			string apath = to_string(result.align.apath);
			printf("\t%d\t%d\t%s\t%d\t%d", 
						result.score,
						result.align.beginPos, apath.c_str(), apath_read_length(result.align.apath), apath_ref_length(result.align.apath));
		}

		printf("\n");


		// save last references, switch buffers if did not use last
		last_refs = refs;		
		if ( !last_refs_used ) {
			linebuf_index = 1 - linebuf_index;
		}
	}

	return 0;
	
}
