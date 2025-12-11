#!/bin/bash

g++ -std=c++17 -g -O0 -Werror \
 -Wno-#pragma-messages \
 -I /Users/drorkessler/Documents/GitHub/homebrew/include \
 \
jbetter.cpp \
 -lhts -lparasail \
 -o jbetter

