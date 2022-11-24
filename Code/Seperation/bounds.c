/**
 * @file   bounds.c
 * @author T. Buckingham
 * @date   Thu May  5 18:05:53 2022
 * 
 * @brief  Functions for finding phoneme boundaries in an audio signal
 * 
 */
#include "bounds.h"

float get_entropy(short* sequence, int inc);

/** 
 * @brief Creates a new array and moves the remaining values to the new array
 * 
 * @param sequence The sequence to procces
 * @param length The length of @param sequence
 * @param shift The amount of data indexed from 0 to remove
 * 
 * @return The new array with the data from 0 -> @param shift
 *
 * Once a boundary has been found and processed this data needs to be removed
 * from the array, this removes this data and returns a new array to save memory
 * and prevent a memory leak.
 */
short* shift_and_reduce(short* sequence, int length, int shift)
{
	short* new = (short*)calloc((length - shift), sizeof(short));
	for(int i = shift; i < length; i++) {
		new[i - shift] = sequence[i];			  
	}
	free(sequence);
	return new;
}

/** 
 * @brief Finds and returns the next boundary in audio signal
 * 
 * @param sequence The audio signal to process
 * @param length The length of @param sequence
 * 
 * @return The next boundary found in the audio sequence
 *
 * The boundary is found using entropy, zero cross and short time energy change
 * thresholds which are defined by the user when running the program.
 */
int next_boundary(short* sequence, int length)
{
	if(length <= 32)
		return -1;
	int inc = 32, gap = 0;
	float change = 0, last_zc = 0, ste_change = 0, last_ste = 0; 
	float entropy = 0, prev_entropy = 0, entropy_change = 0;
	
	for(int i = inc; i < length - inc; i+=inc) {
		entropy = 0;
		float zc = cross_rate(&sequence[i], inc);
		change = abs(last_zc - zc);

		float ste = short_time_energy(&sequence[i], inc, inc);
		if(last_ste == 0) {
			last_ste = 1;
		}
		ste_change = ((ste-last_ste) / fabs(last_ste)) * 100;
		
		
		entropy = get_entropy(&sequence[i], inc); 
		if(prev_entropy == 0) {
			entropy_change = fabs(entropy - prev_entropy);
		} else {
			entropy_change = ((entropy-prev_entropy) / fabs(prev_entropy)) * 100;
		}
		
		
		if((i != 32) && (gap > 6) && 
		   ((change >= glbl_zc_incr && ste_change <= glbl_ste_incr && entropy_change <= glbl_entr_incr) ||
		   (ste_change > glbl_larg_ste_incr && entropy_change > glbl_larg_entr_incr))) {
			return i;
		}
		gap++;
		prev_entropy = entropy;
		last_zc = zc;
		last_ste = ste;
		change = 0;
	}
	return -1;
}

/** 
 * @brief Calculates the entropy value of a given audio chunk
 * 
 * @param sequence The audio chunk to process
 * @param inc The length of @param sequence
 * 
 * @return The entropy value of @param sequence
 */
float get_entropy(short* sequence, int inc)
{
	float entropy = 0;
	for(int m = 0; m < inc; m++) {
		if(sequence[m] == 0) {
			sequence[m] = 1;
		}
		sequence[m] = abs(sequence[m]);
		if(sequence[m] == -32768) {
			sequence[m] = 32767;
		}
		entropy += (float)(sequence[m]) * (log(sequence[m]) / log(2));
	}
	return entropy;
}

