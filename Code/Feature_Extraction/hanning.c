/**
 * @file   hanning.c
 * @author T. Buckingham
 * @date   Thu May  5 17:41:08 2022
 * 
 * @brief  Functions for creating and applying a Hanning window
 *
 */
#include "hanning.h"

/** 
 * @brief Produces a Hanning window of the given length
 * 
 * @param num The length of the Hanning window to produce
 * 
 * @return The Hanning window of length @param num
 */
float* hanning_window(int num)
{
	float* window = (float*)malloc(num * sizeof(float));
	if(window == NULL) 
		printf("Failed to malloc 'window' in hanning\n");
	for(int i = 0; i < num; i++) {
		window[i] = 0.5 * (1 - cos(2*M_PI*(i/num)));
	}

	return window; 
}

/** 
 * @brief Applies a Hanning window to each given chunk of an audio input
 * 
 * @param input The audio signal to be processed
 * @param length The length of @param length
 * @param incr The window width to use
 * 
 * @return The array of Hanning windows
 *
 * The input signal will be blocked into chunks of size @param incr then
 * each chunk will have a Hanning window applied.
 */
float** hanning_chunks(float* input, int length, int incr)
{
	int last = floor((length - (glbl_window_width)) / (glbl_window_width / glbl_interval_div));
	if(last > glbl_frame_limit) {
		last = glbl_frame_limit;
	}
	
	float** chunks = (float**)malloc( (last) * sizeof(float*));
	for(int i = 0; i < (last); i++) {
		chunks[i] = (float*)calloc(incr, sizeof(float));
	}
	int n = 0;
	for(int i = 0; i < (last); i++) {
		float* temp = (float*)malloc(incr * sizeof(float));
		int first = (i) * (incr / glbl_interval_div);
		memcpy(temp, &input[first], incr * sizeof(float));
		float* hanned = hanning(temp, incr);
		memcpy(chunks[n], hanned, sizeof(float) * incr);
		free(hanned);
		n++;
		free(temp);
	}

	return chunks;
}

/** 
 * @brief Applies a Hanning window to each given chunk of an audio input with
 * no overlaping of the windows
 * 
 * @param input The audio signal to be processed
 * @param length The length of @param length
 * @param incr The window width to use
 * 
 * @return The array of Hanning windows
 *
 * The input signal will be blocked into chunks of size @param incr then
 * each chunk will have a Hanning window applied.
 */
float** hanning_chunks_no_overlap(float* input, int length, int incr)
{
	int last = floor(length / incr);
	float** chunks = (float**)malloc( (last) * sizeof(float*));
	for(int i = 0; i < (last); i++) {
		chunks[i] = (float*)calloc(incr,  sizeof(float));
	}
	int n = 0;
	for(int i = 0; i < (last); i++) {
		float* temp = (float*)malloc(incr * sizeof(float));
		int first = (i) * (incr);
		memcpy(temp, &input[first], incr * sizeof(float));
		float* hanned = hanning(temp, incr);
		memcpy(chunks[n], hanned, sizeof(float) * incr);
		free(hanned);
		n++;
		free(temp);
	}

	return chunks;
}

/** 
 * @brief Applies a Hanning function to a given signal
 * 
 * @param array The input signal to apply the Hanning function to
 * @param num The width of the window
 * 
 * @return The result of applying the Hanning function to the given input signal
 */
float* hanning(float* array, int num)
{
	float* results = (float*)calloc(num, sizeof(float));
	float mult = 0;
	for(int i = 0; i < num; i++) {
		mult = 0.5 * (1 - cos(2*M_PI*(i)/(num)));
		results[i] = mult * array[i];
	}
	
	return results;
}
