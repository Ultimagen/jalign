#!/bin/bash -xv

g++ -O3 -Werror \
 -I /manta/manta-1.6.0/src/c++/lib/alignment/ \
 -I /manta/manta-1.6.0/src/c++/lib \
 -DNDEBUG \
 \
 /manta/manta-1.6.0/src/c++/lib/blt_util/blt_exception.cpp \
 /manta/manta-1.6.0/src/c++/lib/blt_util/align_path.cpp \
 /manta/manta-1.6.0/src/c++/lib/blt_util/parse_util.cpp \
 jump_align.cpp \
 -o jump_align

g++ -O3 -Werror \
 -DNDEBUG \
 \
 para_jalign.cpp \
 -o para_jalign \
 -lparasail
 

