#ifndef KNN_H
#define KNN_H

#include <memory.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>

/** 
 * \struct Guess
 * \brief the Guess struct is used to store a guess in KNN for sorting and deciding the result
 * @guess index guessed, can be the phoneme, group or voice
 * @diff the difference result used for sorting
 * @ref_indx the MFCC index used to produce the result - used for storing the amount of correct and incorrect guesses for the MFCC
 */ 
struct Guess {
	int guess;
	double diff;
	int ref_indx;
	struct Phoneme* ref;
};

#include "../Dynamic_Time_Warping/dtw.h"
#include "cluster.h"

int knn_mfccs(float** test, int test_length, int k, char* ph);
int knn_mfccs_size(float** total_test, int test_length, int k, char* ph);
int k_means(float* test, int test_length, int k);
int knn_mfccs_size_noref(float** test, int test_length, int k, char* ph);
int knn_mfccs_voice_time(float* test, int test_length, int k, char* ph);
int knn_mfccs_group_time(float* test, int test_length, int k, char* ph);

#endif
