#ifndef DTW_H
#define DTW_H

#include "../Misc/includes.h"

struct Ph_index {
	int i;
	char group[5];
	char* name;
};

struct Phoneme {
	struct Ph_index* index;
	double score;
	float** mfcc;
	int* size;
	int* count;
	int size_count;
	int use_count;
} ph;

long double dtw_frame_result(float* signal, int signal_length, struct Phoneme* phoneme, int p);

#endif
