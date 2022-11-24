/**
 * @file   mel.c
 * @author T. Buckingham
 * @date   Thu May  5 14:56:25 2022
 * 
 * @brief  Functions for generating a Mel filter bank and applying it to an FFT'ed sequence.
 *
 */
#include "mel.h"

/** 
 * @brief Generates a Mel scale filter bank
 * 
 * @param width The width of the FFT'ed input sequence
 * @param banks The number of banks to be generated
 * 
 * @return The Mel filter bank of size @param width by @param banks
 *
 * Generates a Mel filter bank with minimum and maximum values of 0 and 4000 respecitively.
 * @var glbl_nfft is used when creating the filter bank.
 */
float** mel_fb(int width, int banks)
{
	float minmel   = 2595.0*log10(1.0 + ( 0 /700.0));
	float maxmel   = 2595.0*log10(1.0 + ( 4000 /700.0));     
	float binwidth = (maxmel - minmel) / (banks + 1);
	// float L = (maxmel - minmel) / (banks * (1 - 0.4) + 0.4);
	float start = minmel; 
	float* mel = calloc(banks + 2, sizeof(float));
	float* h = calloc(banks + 2, sizeof(float)); 
	float* f = calloc(banks + 2, sizeof(float));
	
	float** H = (float**)malloc(width * sizeof(float*));
	if(H == NULL)
		printf("Failed to malloc H in mel_fb\n");
	for(int i = 0; i < width; i++) {
		H[i] = (float*)calloc(banks + 2, sizeof(float));
		if(H[i] == NULL)
			printf("Failed to calloc H[i] in mel_fb\n");
		for(int n = 0; n < banks + 2; n++) {
			H[i][n] = 0;
		}
	}
	
	for(int n = 0; n < banks + 2; n++) {
		mel[n] = start; 
		h[n] = 700.0*(exp(mel[n]/1125.0) - 1);
		f[n] = floor((glbl_nfft + 1) * h[n] / 16000); // 32, 64, 128, 256, 512 (orignally using 128)
		// float binwidth = (L / 2) + n * (maxmel - minmel - L) / (banks - 1) + minmel;  
		start = start + binwidth;
	}
	
	/* /\* TODO test with m = 2, k = 1, [m], [m - 1], [k], [k - 1], etc.*\/ */
	/* for(float m = 2; m < banks + 2; m++) { */
	/* 	for(float k = 0; k < width; k++) { */
	/* 		if (k >= f[(int)m - 1] && k <= f[(int)m] && (f[(int)m] - f[(int)m - 1] != 0)) { */
	/* 			H[(int)k][(int)m - 2] = (k - f[(int)m - 1]) / (f[(int)m] - f[(int)m - 1]); */
	/* 		} else if (k >= f[(int)m] && k <= f[(int)m + 1] && (f[(int)m + 1] - f[(int)m] != 0)) { */
	/* 			H[(int)k][(int)m - 2] = (f[(int)m + 1] - k) / (f[(int)m + 1] - f[(int)m]); */
	/* 		} else { */
	/* 			H[(int)k][(int)m - 2] = 0; */
	/* 		} */
	/* 		if(fpclassify(H[(int)k][(int)m - 2]) == FP_INFINITE || fpclassify(H[(int)k][(int)m - 2]) == FP_NAN) { */
	/* 			H[(int)k][(int)m - 2] = 0; */
	/* 		} */
	/* 	} */
		
	/* } */
	
	for(float m = 1; m < banks + 1; m++) {
		for(float k = 0; k < width; k++) {
			if (k >= f[(int)m - 1] && k <= f[(int)m] && (f[(int)m] - f[(int)m - 1] != 0)) {
				H[(int)k][(int)m - 1] = (k - f[(int)m - 1]) / (f[(int)m] - f[(int)m - 1]);
			} else if (k >= f[(int)m] && k <= f[(int)m + 1] && (f[(int)m + 1] - f[(int)m] != 0)) {
				H[(int)k][(int)m - 1] = (f[(int)m + 1] - k) / (f[(int)m + 1] - f[(int)m]);
			} else {
				H[(int)k][(int)m - 1] = 0;
			}
			if(fpclassify(H[(int)k][(int)m - 1]) == FP_INFINITE || fpclassify(H[(int)k][(int)m - 1]) == FP_NAN) {
				H[(int)k][(int)m - 1]= 0;
			}
		}
	}

	/* for(int k = 0; k < width; k++) { */
	/* 	for(int m = 0; m < banks + 1; m++) { */
	/* 		printf("%f, ", H[k][m]); */
	/* 	} */
	/* 	printf("; \n"); */
	/* } */
	/* exit(0); */
	free(mel);
	free(h);
	free(f);
	
	return H;
}

/** 
 * @brief Applies a generated Mel filter bank to the input array
 * 
 * @param array The FFT'ed sequence the filter bank is to be applied to
 * @param width The length of @param array
 * @param banks The number of desired filter banks
 * 
 * @return The result of applying the Mel filter bank to the input FFT'ed sequence.
 */
float* mel(float* array, int width, int banks)
{
	float* result = (float*)calloc(banks, sizeof(float));
	float** melfb = mel_fb(width, banks);
	float** applied = calloc(width, sizeof(float*));
	if(result == NULL || melfb == NULL || applied == NULL)
		printf("Malloc/Calloc error in mel\n");
	for(int i = 0; i < width; i++) {
		applied[i] = calloc(banks, sizeof(float));
		if(applied[i] == NULL)
			printf("Failed to calloc applied[i]\n");
		for(int j = 0; j < banks; j++) {
			applied[i][j] = 0;
		}
	}
	
	for(int i = 0; i < width; i++) {
		for(int j = 0; j < banks; j++) {
			applied[i][j] = array[i] * melfb[i][j];
			if(fpclassify(melfb[i][j]) == FP_INFINITE || fpclassify(melfb[i][j]) == FP_NAN) {
				printf("MelFb result contained NAN of INF :: %.3f || width := %d || i := %d || j := %d\n", melfb[i][j],  width, i, j);
	 		}
			if(fpclassify(array[i]) == FP_INFINITE || fpclassify(array[i]) == FP_NAN) {
				printf("Array contained NAN of INF :: %.3f || width := %d\n", melfb[i][j],  width);
			}
			if(applied[i][j] == 0) {
				applied[i][j] = FLT_EPSILON;
			}
		}
	}
	
	for(int i = 0; i < width; i++) {
		for(int j = 0; j < banks; j++) {
			result[j] += applied[i][j];
		}		
	}
	// maybe change
	/* for(int i = 0; i < banks; i++) { */
	/* 	result[i] = result[i] / banks; */
	/* } */
	for(int i = 0; i < width; i++) {
		free(melfb[i]);
	}
	free(melfb);
	for(int i = 0; i < width; i++) {
		free(applied[i]);
	}
	free(applied);
	
	return result;
}

	     
