#ifndef CROSS_RATE_H
#define CROSS_RATE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#include "../Feature_Extraction/fft.h"
#include "../Feature_Extraction/hanning.h"
#include "../Dynamic_Time_Warping/dtw.h"

int cross_rate(short* signal, int signal_length);
float f_cross_rate(float* signal, int signal_length);
float stavg_cross_rate(float* signal, int signal_length);
float stavg_cross_rate_no_overlap(float* signal, int signal_length);
float favg_cross_rate(short* signal, int signal_length);

#endif

