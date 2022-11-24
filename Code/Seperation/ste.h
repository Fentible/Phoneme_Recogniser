#ifndef STE_H
#define STE_H

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <float.h>

#include "../Feature_Extraction/hanning.h"

float short_time_energy(short* signal, int signal_length, int window_length);
float f_short_time_energy(float* signal, int signal_length, int window_length);

#endif
