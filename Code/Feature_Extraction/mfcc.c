/**
 * @file   mfcc.c
 * @author T. Buckingham
 * @date   Thu May  5 16:06:24 2022
 * 
 * @brief  Functions for creating MFCCs and additional MFCC features, and normalising MFCCs.
 * 
 */
#include "mfcc.h"

float min_mfcc = FLT_MAX;  /* The minimum value found in MFCCs during training; used for normalisation */
float max_mfcc = FLT_MIN;  /* The maximum value found in MFCCs during training; used for normalisation */

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

/** 
 * @brief Calculates the log entropy of a given signal.
 * 
 * @param sequence The signal to be processed.
 * @param inc The length of @param sequence.
 * 
 * @return The log entropy value.
 */
float log_entropy(float* sequence, int inc)
{
	long double entropy = 0;
	for(int m = 0; m < inc; m++) {
		if(sequence[m] == 0) {
			sequence[m] = 1;
		}
		sequence[m] = fabs(sequence[m]);
		entropy += (sequence[m]) * (log(sequence[m]) / log(2));
	}
	if(entropy == 0.0)
		return 0;
	return (float)log(fabs(entropy));
}

/** 
 * @brief Calculates the log energy of a given signal.
 * 
 * @param sequence The signal to be processed.
 * @param inc The length of @param sequence.
 * 
 * @return The log energy value.
 */
float log_energy(float* chunk, int length)
{
	long double sum = 0;
	for(int i = 0; i < length; i++) {
		sum += chunk[i] * chunk[i];
	}
	if(sum == 0)
		return 0;
	else 
		return log10(sum);
}

/** 
 * @brief Calculates the short time energy of a given signal.
 * 
 * @param sequence The signal to be processed.
 * @param inc The length of @param sequence.
 * 
 * @return The short time energy value.
 */
float st_energy(float* chunk, int length)
{
	long double sum = 0;
	for(int i = 0; i < length; i++) {
		sum += chunk[i] * chunk[i];
	}
	return sum;
}

/** 
 * @brief Calculates the kurtosis of a given signal.
 * 
 * @param sequence The signal to be processed.
 * @param inc The length of @param sequence.
 * 
 * @return The log kurtosis value.
 */
float kurtosis(float* chunk, int length)
{
	float result = 0, mean = 0, top = 0, bot = 0;
	float N = 1.0/length;
	for(int i = 0; i < length; i++) {
		mean += chunk[i];
	}
	mean /= length;
	for(int i = 0; i < length; i++) {
		if(fabs(chunk[i] - mean) <= FLT_EPSILON) {
			top += 0;
			bot += 0;
		} else {
			top += pow( (chunk[i] - mean), 4);
			bot += pow( (chunk[i] - mean), 2);
		}
	}
	if(N == 0 || bot == 0) {
		return 0;
	}
	result = (N * top) / pow((N * bot), 2);
	return result;
}

/** 
 * @brief Normalises an MFCC using @var min_mfcc and @var max_mfcc
 * 
 * @param mfcc The MFCC sequence to be normalised.
 * @param length The length of @param mfcc
 */
void normalise_mfcc(float* mfcc, int length)
{
	for(int i = 0; i < length; i++) {
		mfcc[i] = (mfcc[i] - min_mfcc) / (max_mfcc - min_mfcc) * 5;
	}
	return;
}

/** 
 * @brief Updates @var min_mfcc and @var max_mfcc if new minimum and maximum values are found.
 * 
 * @param mfcc The MFCC sequence to be processed.
 * @param length The length of @param mfcc
 */
void update_mfcc_norm(float* mfcc, int length)
{
	for(int i = 0; i < length; i++) {
		if(mfcc[i] > max_mfcc)
			max_mfcc = mfcc[i];
		if(mfcc[i] < min_mfcc)
			min_mfcc = mfcc[i];
	}
	return;
}

/** 
 * @brief The control function that generates an MFCC from the given input audio sequence.
 * 
 * @param sequence The audio sequence to be processed
 * @param width The length of @param sequence
 * @param incr The window width used for Hanning
 * @param banks The number of Mel filter banks to use
 * @param paa The Piece-wise aggregation divisor
 * 
 * @return The generated MFCC.
 */
float* mfcc(float* sequence, int width, int incr, int banks, int paa)
{

	if(paa == 0) {
		sequence = f_paa(sequence, width);
		width = floor(width / glbl_paa);
	}
	/* if width < width * 2 || width * 1.5? then use non overlapping */
	int amount = floor((width - (incr)) / (incr / glbl_interval_div));
        if(amount > glbl_frame_limit) {
		amount = glbl_frame_limit;
	}
	/* all points in window applied to all points of the window function */
	float** chunks = hanning_chunks(sequence, width, incr);
	if(amount <= 0) {
		free(chunks);
		chunks = NULL;
	}
	if(chunks == NULL) {
		// printf("Hanning chunks with overlaps failed trying no overlap :: \n");
		chunks = hanning_chunks_no_overlap(sequence, width, incr);
		amount = floor(width / incr);
		if(chunks == NULL) {
			printf("Oops, no overlap failed as well, sequence is too short or window width is too large\nNo MFCC will be created or tested :: width : %d || incr :: %d\n", width, incr);
			return NULL;
		}
	}

	
	float** mags = fft_chunks(chunks, amount, incr);
	
	for(int i = 0; i < amount; i++) { 
		for (int low = 0, high = (incr / 2) - 1; low < high; low++, high--) {
			float temp = mags[i][low];
			mags[i][low] = mags[i][high];
			mags[i][high] = temp;
		}
	}
	
	float** applied_mels = calloc(amount, sizeof(float*));
	if(chunks == NULL || mags == NULL || applied_mels == NULL) {
		printf("Malloc/Calloc error in mfcc\n");
		return NULL;
	}
		
	for(int i = 0; i < amount; i++) {
		applied_mels[i] = (float*)calloc(banks,  sizeof(float));
		if(applied_mels[i] == NULL) {
			printf("Failed to calloc applied_mels[i]\n");
			return NULL;
		}
		float* temp = mel(mags[i], (incr / 2), banks);
		memcpy(applied_mels[i], temp, sizeof(float) * banks);
		free(temp);
	}
	for(int i = 0; i < amount; i++) {
		for(int j = 0; j < banks; j++) {
			/* TODO :: Conditional jump or move depends on uninitialised value(s) */
			if(applied_mels[i][j] <= 0) {
				applied_mels[i][j] = FLT_EPSILON;
			}
			applied_mels[i][j] = log10(applied_mels[i][j]);	
		}
	}
	for(int i = 0; i < amount; i++) {
		float* temp = dct(applied_mels[i], banks);
		for(int j = 0; j < banks; j++) {
			applied_mels[i][j] = temp[j];
			// memcpy(applied_mels[i], temp, sizeof(float) * banks);
		}
		free(temp);
	}
	int trunc = floor(banks * glbl_test_trunc);
	// int trunc = 1;
	float* result = (float*)malloc((amount * trunc) * sizeof(float));
	
	int last = amount;

	int n = 0;
	int m = 0;
	for(int i = 0; i < amount; i++) {
		for(int j = 0; j < trunc; j++) {
			if(fpclassify(applied_mels[i][j]) == FP_INFINITE || fpclassify(applied_mels[i][j]) == FP_NAN) {
				result[n] = 0; // FLT_EPSILON;
				printf("MFCC result contained NAN of INF :: %.3f || trunc := %d || amount := %d || i := %d || j := %d\n", applied_mels[i][j], trunc, amount, i, j);
			} else {
				result[n] = applied_mels[i][j];
			}
			if(LOG_E) {
				if(m == (trunc - 1)) {
					result[n] = f_cross_rate(chunks[i], glbl_window_width);
				}
				m++;
				if(m >= trunc) {
					m = 0;
				}
				
			}
			n++;
		}
		
		free(applied_mels[i]);
	}
	free(applied_mels);

	for(int i = 0; i < (last); i++) {
		free(chunks[i]);
	}
	free(chunks);

	for(int i = 0; i < amount; i++) {
		free(mags[i]);
	}
	free(mags);
	free(sequence);	
	// mfcc_window(result, trunc, amount);
	return result;
}

