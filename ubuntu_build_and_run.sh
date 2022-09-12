#!/bin/bash

OPTARG=0
CONDITION=${1:-0}
if [ $# != 0 ]
then
OPTARG=3
fi

rm executable.exe
g++ -O$OPTARG -g -std=c++11 -Wall -Wextra -Wno-misleading-indentation -Wno-char-subscripts -Wconversion -Wno-write-strings -Wno-missing-field-initializers -o executable.exe hw.cpp -lglfw -lGL
./executable.exe
