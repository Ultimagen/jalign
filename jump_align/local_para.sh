#!/bin/bash

g++ -std=c++11 -g -O3 -Werror \
 -Wno-#pragma-messages \
 -I /Users/drorkessler/Documents/GitHub/homebrew/include \
 \
para_jalign.cpp \
 -lparasail \
 -o para_jalign

