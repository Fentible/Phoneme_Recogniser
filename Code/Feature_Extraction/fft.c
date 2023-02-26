/**
 * @file   fft.c
 * @author T. Buckingham
 * @date   Thu May  5 14:48:03 2022
 * 
 * @brief Fast Fourier transform and Discrete Cosine tranform functions 
 * 
 */
#include "fft.h"

/** 
 * @brief DCT function used during the MFCC process
 * 
 * @param array The MFCC to be processed
 * @param width The length of @param array
 * 
 * @return 
 */
float* dct(float* array, int width)
{
	float* result = (float*)malloc(width * sizeof(float));
	if(result == NULL)
		printf("Failed to malloc 'result' in DCT\n");
	
	for(int i = 0; i < width; i++) {
		float sum = 0;
		for(int n = 0; n < width; n++) {
			sum += array[n] * cos(M_PI / (width) * (n + 0.5) * i);
			
		}
		result[i] = sum;
	}
	return result;
}

/** 
 * @brief Produces an 2 dimensional array of FFT magnitudes
 * 
 * @param chunks The Hanning window sequences from an audio signal
 * @param n The number of Hanning windows to be processed
 * @param length The length of each Hanning window
 * 
 * @return A 2D array of FFT arrays
 *
 * Each 'chunk' passed will have an FFT applied and the magnitude extracted, the same amount of arrays will be returned
 * however, they will all be of length (@param length / 2)
 */
float** fft_chunks(float** chunks, int n, int length)
{
	float complex* fft_buff = calloc(length, sizeof(float complex));
	float complex* result;
	float** mag = (float**)calloc(n, sizeof(float*));
	if(mag == NULL)
		printf("Failed to malloc 'mag' in fft\n");
	for(int i = 0; i < n; i++) {
		mag[i] = (float*)calloc((length / 2), sizeof(float));
		if(mag[i] == NULL)
			printf("Failed to malloc 'mag[i]' in fft\n");
	}
	for(int m = 0; m < n; m++) {
		for(int i = 0; i < length; i++) {
			fft_buff[i] = chunks[m][i];
		}
		result = fft(fft_buff, length);
		for(int i = 0; i < (length / 2); i++) {
			/* TODO :: Conditional jump or move depends on uninitialised value(s) */
			mag[m][i] = cabsf(result[i]); //sqrt( (real[i] * real[i]) + (imag[i] * imag[i]) );
		}
		free(result);
	}

	free(fft_buff);
	return mag;
}

/** 
 * @brief Performs an FFT on a given array
 * 
 * @param chunk The complex array to be processed
 * @param length The length of @param chunk
 * 
 * @return The FFT in the form of a complex number, no magnitude taken
 *
 * Before processing the input must be in complex form and the return will be also be in complex form and so
 * the magnitude will need to be taken for further MFCC processing.
 */
float complex* fft(float complex* chunk, int length)
{
	int n = length; 
	if (n == 1)
		return chunk;

	double complex w = cexp( (2 * M_PI * I) / n);
	float complex* even = calloc(length / 2, sizeof(float complex));
	float complex* odd = calloc(length / 2, sizeof(float complex));
	if(odd == NULL || even == NULL) {
		printf("Failed to malloc 'odd || even' in FFT\n");
	}
	int e = 0, o = 0;
	for(int i = 0; i < length; i++) {
		if(i % 2 == 0) {
			even[e] = chunk[i];
			e++;
		} else {
			odd[o] = chunk[i];
			o++;
		}
	}
	float complex* y_even = fft(even, length / 2);
	float complex* y_odd = fft(odd, length / 2);
	float complex* y = (float complex*)malloc(length * sizeof(float complex));
	if(y == NULL) {
		printf("Failed to malloc 'y' in FFT\n");
	}
	for(int j = 0; j < length / 2; j++) {
		y[j] = y_even[j] + (cpow(w, j) * y_odd[j]);
		y[j + n / 2] = y_even[j] - (cpow(w, j) * y_odd[j]);
	}
	
	if(even != y_even)
		free(y_even);
	if(odd != y_odd)
		free(y_odd);
	free(even);
	free(odd);
	return y;
}
