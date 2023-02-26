#ifndef DELTA_H
#define DELTA_H

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "../Clustering/cluster.h"
#include "../Dynamic_Time_Warping/dtw.h"

float* delta(float* mfcc, int size);
float* delta_delta(float* mfcc_delta, int size);
void update_deltas_norm(float* delta, int length);
void update_delta_deltas_norm(float* delta, int length);
void normalise_delta_delta(float* delta, int length);
void normalise_delta(float* delta, int length);

#endif
