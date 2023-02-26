#ifndef MEL_H
#define MEL_H

#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "../Clustering/cluster.h"
#include "../Dynamic_Time_Warping/dtw.h"

float** mel_fb(int width, int banks);
float* mel(float* array, int width, int banks);

#endif
