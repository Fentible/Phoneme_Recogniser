#ifndef MFCC_H
#define MFCC_H

#include "mel.h"
#include "paa.h"
#include "hanning.h"
#include "fft.h"

#include "../Clustering/cluster.h"

float* mfcc(float* sequence, int size, int window, int banks, int paa);
float log_energy(float* chunk, int length);
void update_mfcc_norm(float* mfcc, int length);
void normalise_mfcc(float* mfcc, int length);
float kurtosis(float* chunk, int length);
float log_entropy(float* sequence, int inc);

#endif
