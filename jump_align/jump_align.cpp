#include <stdio.h>
#include <string.h>
#include <string>
#include "GlobalAligner.hpp"
#include "GlobalJumpAligner.hpp"
#include "logging.h"

// some hard limits
#define SEQ_MAX_LEN	5000
#define REF_MAX_LEN 5000

// calculated limits
#define LINE_MAX_LEN (SEQ_MAX_LEN + 1 + REF_MAX_LEN + 1 + REF_MAX_LEN + 1 + 1)

// line buffers used for reading input
// double buffering used to optimize the (common?) case where the reference repeats from the previous line
static char linebuf[2][LINE_MAX_LEN];
static int linebuf_index = 0;

// other local vars
static const char* valid_bases = "ACGTNacgtn";

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
	if ( argc != 7 ) {
		fprintf(stderr, "wrong score arguments:\n");
		fprintf(stderr, "usage: %s <match> <mismatch> <open> <extend> <offedge> <jump>\n", argv[0]);
		exit(-1);
	}

	// configure aligners
	AlignmentScores<int> scores(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), false);
	GlobalJumpAligner<int> jump_aligner(scores, atoi(argv[6]));
	GlobalAligner<int> aligner(scores);

	// loop on reading lines and process
	int lineno = 0;
	struct refs_t last_refs;
	printf("score\tjumpInsertSize\tjumpRange\tbeginPos1\tapath1\tbeginPos2\tapath2\tscore1\tbegin1\tcigar1\tscore2\tbegin2]tcigar2\n");
	while ( fgets(linebuf[linebuf_index], sizeof(linebuf[0]), stdin) ) {
		lineno++;
		char* line = linebuf[linebuf_index];

		// ignore empty or comment lines
		if ( !line[0] || line[0] == '\n' || line[0] == '#' ) {
			continue;
		}

		// ignore lines not starting with a valid base
		if ( !strchr(valid_bases, line[0]) ) {
			continue;
		}

		// parse the line, get sequence
		char* seq = strtok(line, "\t");
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

		// jump align
		JumpAlignmentResult<int> result;
		jumpAlign<GlobalJumpAligner<int>,const char*,int>(
				jump_aligner, seq, seq + seq_len, refs.ptr[0], refs.ptr[0] + refs.len[0], refs.ptr[1], refs.ptr[1] + refs.len[1], result);
		string apath1 = to_string(result.align1.apath);
		string apath2 = to_string(result.align2.apath);
		printf("%d\t%d\t%d\t%d\t%s\t%d\t%s", 
					result.score, result.jumpInsertSize, result.jumpRange,
					result.align1.beginPos, apath1.c_str(),
					result.align2.beginPos, apath2.c_str()
					);

		// simple align
		for ( int i = 0 ; i < 2 ; i++ ) {
			AlignmentResult<int> result;
			aligner.align(seq, seq + seq_len, refs.ptr[i], refs.ptr[i] + refs.len[i], result);
			string apath = to_string(result.align.apath);
			printf("\t%d\t%d\t%s", 
						result.score,
						result.align.beginPos, apath.c_str()
						);
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
