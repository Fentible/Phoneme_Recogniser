#ifndef REALLOC_H
#define REALLOC_H

#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

short* s_realloc(short* input, int length);
float* f_realloc(float* input, int length);
int* i_realloc(int* input, int length);

#endif
