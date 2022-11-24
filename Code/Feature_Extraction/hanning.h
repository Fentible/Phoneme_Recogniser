#ifndef HANNING_H
#define HANNING_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "../Training/train.h"
#include "../Seperation/cross_rate.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

float* hanning(float* array, int num);
float* hanning_window(int num);
float** hanning_chunks_no_overlap(float* input, int length, int incr);
float** hanning_chunks(float* input, int length, int incr);

#endif
