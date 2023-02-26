#ifndef MFCC_H
#define MFCC_H

#include "../Misc/includes.h"

int frame_amount(int signal_length);
int mfcc_size(int signal_length);
float* mfcc_quick(float* sequence, int width, int incr, int banks, int paa);

#endif
