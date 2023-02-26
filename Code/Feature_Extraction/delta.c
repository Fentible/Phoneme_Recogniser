/**
 * @file   delta.c
 * @author T. Buckingham
 * @date   Wed May  4 23:15:48 2022
 * 
 * @brief  Functions for creating delta and detla-delta coefficients of MFCCs
 * 
 * 
 */
#include "delta.h"

float min_delta = FLT_MAX;         /* Minimum delta value used for normalisation */
float max_delta = FLT_MIN;         /* Maximum delta value used for normalisation */

float min_delta_delta = FLT_MAX;   /* Minimum delta-delta value used for normalisation */
float max_delta_delta = FLT_MIN;   /* Minimum delta-delta value used for normalisation */

/** 
 * @brief Normalises a delta coefficient array
 * 
 * @param delta The delta array to be normalised
 * @param length The length of @delta
 */
void normalise_delta(float* delta, int length)
{
	for(int i = 0; i < length; i++) {
		delta[i] = (delta[i] - min_delta) / (max_delta - min_delta) * 3;
	}
	return;
}

/** 
 * @brief Normalises a detal-delta coefficient array
 * 
 * @param delta The delta-delta coefficients array to be normalised
 * @param length The length of @param delta
 */
void normalise_delta_delta(float* delta, int length)
{
	for(int i = 0; i < length; i++) {
		delta[i] = (delta[i] - min_delta) / (max_delta - min_delta);
	}
	return;
}

/** 
 * @brief Update the minimum and maximum values of the delta coefficients 
 * 
 * @param delta The delta sequence to process
 * @param length The length of @param delta
 */
void update_deltas_norm(float* delta, int length)
{
	for(int i = 0; i < length; i++) {
		if(delta[i] > max_delta)
			max_delta = delta[i];
		if(delta[i] < min_delta)
			min_delta = delta[i];
	}
	return;
}

/** 
 * @brief Update the minimum and maximum values of the delta-delta coefficients 
 * 
 * @param delta The delta-delta sequence to process
 * @param length The length of @param delta-delta
 */
void update_delta_deltas_norm(float* delta_delta, int length)
{
	for(int i = 0; i < length; i++) {
		if(delta_delta[i] > max_delta_delta)
			max_delta = delta_delta[i];
		if(delta_delta[i] < min_delta)
			min_delta = delta_delta[i];
	}
	return;
}

/** 
 * @brief Calculate the delta coefficients of from an MFCC
 * 
 * @param mfcc The MFCC to process
 * @param size The length of @param MFCC
 * 
 * @return The delta coefficient array
 */
float* delta(float* mfcc, int size)
{

	float* delta = (float*)calloc(size, sizeof(float));
	int trunc = floor((glbl_banks) * glbl_test_trunc);
	// k1 is forward and k2 is backwards
	int signal_length = floor(size / trunc);
	int k1 = 0, k2 = 0;
	for(int i = 0; i < signal_length; i++) {
		k1 = fmin(abs(i - (signal_length - 1)), 2); 
		k2 = fmin(abs(0 - i), 2);
		for(int m = 0; m < trunc; m++) {
			delta[(i * trunc) + m] = mfcc[((i + k1) * trunc) + m] - mfcc[((i - k2) * trunc) + m];
			if(fpclassify(delta[(i * trunc) + m]) == FP_INFINITE || fpclassify(delta[(i * trunc) + m]) == FP_NAN) {
				perror("SIGPE\n");
				raise(SIGFPE);
			}
		}
	}
	return delta;
}

/** 
 * @brief Calculate the delta-delta coefficients of from a delta coefficient array
 * 
 * @param mfcc The delta coefficient array to process
 * @param size The length of @param mfcc_delta
 * 
 * @return The delta coefficient array
 */
float* delta_delta(float* mfcc_delta, int size)
{

	float* delta_delta = (float*)calloc(size, sizeof(float));
	int trunc = floor((glbl_banks) * glbl_test_trunc);
	// k1 is forward and k2 is backwards
	int signal_length = floor(size / trunc);
	int k1 = 0, k2 = 0;
	for(int i = 0; i < signal_length; i++) {
		k1 = fmin(abs(i - (signal_length - 1)), 2); 
		k2 = fmin(abs(0 - i), 2);
		for(int m = 0; m < trunc; m++) {
			delta_delta[(i * trunc) + m] = mfcc_delta[((i + k1) * trunc) + m] - mfcc_delta[((i - k2) * trunc) + m];
			if(fpclassify(delta_delta[(i * trunc) + m]) == FP_INFINITE || fpclassify(delta_delta[(i * trunc) + m]) == FP_NAN) {
				perror("SIGPE\n");
				raise(SIGFPE);
			}
		}
	}
	return delta_delta;
}
