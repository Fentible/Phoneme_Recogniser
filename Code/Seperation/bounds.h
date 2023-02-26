#ifndef BOUNDS_H
#define BOUNDS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <float.h>
#include "cross_rate.h"
#include "ste.h"
#include "../Feature_Extraction/fft.h"

int next_boundary(short* sequence, int length);
short* shift_and_reduce(short* sequence, int length, int shift);
	
#endif
