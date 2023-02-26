/**
 * @file   realloc.c
 * @author T. Buckingham
 * @date   Fri Apr 29 01:07:49 2022
 * 
 * @brief  Safe reallocation functions for arrays
 * 
 * A range of reallocation functions which contain checks before
 * return the new array. Contains a int, short and float variant.
 * 
 */
#include "realloc.h"

/** 
 * @brief Reallocates @input to size @length
 * 
 * @param input The array to reallocate
 * @param length The resulting size of the new array
 * 
 * @return An array of size @length
 */
short* s_realloc(short* input, int length) {

	short* new = realloc(input, sizeof(short) * length);

	if (new != NULL) {
		return new;
	} else {
		fprintf(stderr, "Memory error, trying to allocate %d length which is %d bytes.\n", (int)length, (int)(length*sizeof(short)));
		free(input);
		exit(-1);
	}
	
}

float* f_realloc(float* input, int length) {

	float* new = realloc(input, sizeof(float) * length);
	
	if (new != NULL) {
		return new;
	} else {
		fprintf(stderr, "Memory error, trying to allocate %d length which is %d bytes.\n", (int)length, (int)(length*sizeof(float)));
		free(input);
		exit(-1);
	}
	
}

int* i_realloc(int* input, int length) {

	int* new = realloc(input, sizeof(int) * length);

	if (new != NULL) {
		return new;
	} else {
		fprintf(stderr, "Memory error, trying to allocate %d length which is %d bytes.\n", (int)length, (int)(length*sizeof(short)));
		free(input);
		exit(-1);
	}
	
}
