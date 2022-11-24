/**
  * @file   paa.c
  * @author T. Buckingham
  * @date   Thu May  5 15:50:21 2022
  * 
  * @brief  Functions to apply Piece-wise aggregation to an input signal
  * 
  * 
  */
#include "paa.h"

/** 
 * @brief Applies PAA to the input sequence
 * 
 * @param sequence The input signal to be processed
 * @param length The length of @param input
 * 
 * @return The new processed signal of size @var length / @var glbl_paa
 *
 * The input signal is divided by the PAA amount set by the user, with each grouping being averaged into one value
 */
short* paa(short* sequence, int length) {

	int n_length = floor(length/glbl_paa);
	
	short* temp = malloc(sizeof(short) * length);
	for(int i = 0; i < length; i++) {
		memcpy(&temp[i], &sequence[i], sizeof(short));
	}
	short* new = s_realloc(sequence, n_length);

	if (new == NULL) {
		fprintf(stderr, "Memory error, trying to allocate %d length which is %d bytes.\n", (int)length, (int)(length*sizeof(short)));
		free(sequence);
		exit(-1);
	}
	int sum = 0, num = 0;
	for(int i = 0; i < glbl_paa; i++) {
		sum = sum + temp[i];
		if(i % glbl_paa == 0) {
			num = sum / glbl_paa;
			new[i / glbl_paa] = num;
			num = 0;
			sum = 0;
		}
	}
	free(temp);
	return new;
	
}

/** 
 * @brief Applies PAA to the input sequence
 * 
 * @param sequence The input signal to be processed
 * @param length The length of @param input
 * 
 * @return The new processed signal of size @var length / @var glbl_paa
 *
 * The input signal is divided by the PAA amount set by the user, with each grouping being averaged into one value.
 * The same as @fn paa() except this version is for floats.
 */
float* f_paa(float* sequence, int length) {

	int n_length = floor(length/glbl_paa);
	
	float* temp = malloc(sizeof(float) * length);
	for(int i = 0; i < length; i++) {
		memcpy(&temp[i], &sequence[i], sizeof(float));
	}
	float* new = f_realloc(sequence, n_length);
	
	if (new == NULL) {
		fprintf(stderr, "Memory error, trying to allocate %d length which is %d bytes.\n", (int)length, (int)(length*sizeof(float)));
		free(sequence);
		exit(-1);
	}
	int sum = 0, num = 0;
	for(int i = 0; i < glbl_paa; i++) {
		sum = sum + temp[i];
		if(i % glbl_paa == 0) {
			num = sum / glbl_paa;
			new[i / glbl_paa] = num;
			num = 0;
			sum = 0;
		}
	}
	free(temp);
	return new;
	
}
