#ifndef FFT_H
#define FFT_H

#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>
#include <complex.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif


float** fft_chunks(float** chunks,int n, int length);
float complex* fft(float complex* chunk, int length);
float* dct(float* array, int width);

#endif
