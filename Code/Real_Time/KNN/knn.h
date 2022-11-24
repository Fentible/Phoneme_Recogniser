#ifndef KNN_H
#define KNN_H

#include "../Misc/includes.h"

struct Guess {
	int guess;
	long double diff;
	int ref_indx;
	struct Phoneme* ref;
};

int knn_mfccs_size(float* test, int test_length, int k);

#endif
