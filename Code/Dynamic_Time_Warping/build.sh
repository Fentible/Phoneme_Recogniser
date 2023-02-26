#!/bin/bash

    eval "gcc -O3 -g -std=c11 ./dtw.c ../Training/train.c ../Misc/realloc.c ../Testing/test.c ../Seperation/cross_rate.c  ../Seperation/ste.c ../Feature_Extraction/*.c ../Seperation/bounds.c ../Clustering/cluster.c ../Clustering/knn.c -D_XOPEN_SOURCE=600 -pthread -onan -o dtw.exe -lm;"

