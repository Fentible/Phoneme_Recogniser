#!/bin/bash

    eval "gcc -g -Wall -Werror -pedantic -std=c11 ./Misc/*.c ./DTW/*.c ./MFCCs/*.c ./KNN/*.c ./Boundary/*.c ./Feature/*.c ./Test/*.c ./*.c -D_XOPEN_SOURCE=600 -pthread -onan -o rte.exe -lm;"

