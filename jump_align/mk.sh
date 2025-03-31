#!/bin/bash -xv

g++ -O3 -Werror \
 -I ../manta/src/c++/lib/alignment/ \
 -I ../manta/src/c++/lib \
 -DNDEBUG \
 \
 ../manta/src/c++/lib/blt_util/blt_exception.cpp \
 ../manta/src/c++/lib/blt_util/align_path.cpp \
 ../manta/src/c++/lib/blt_util/parse_util.cpp \
 jump_align.cpp \
 -o jump_align

