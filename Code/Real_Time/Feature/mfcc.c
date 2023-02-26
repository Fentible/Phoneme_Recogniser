#include "mfcc.h"

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

int mfcc_size(int signal_length)
{

	int width = floor(signal_length / glbl_paa);
	int amount = floor((width - (glbl_window_width)) / (glbl_window_width / glbl_interval_div)) + 1;
	if(amount <= 0) {
		amount = floor(width / glbl_window_width);
	}
	if(amount > 15) {
		amount = 15;
	}
	int new_size = amount * floor((glbl_banks) * glbl_test_trunc);
	return new_size;
	
}

int frame_amount(int signal_length)
{

	int width = floor(signal_length / glbl_paa);
	int amount = floor((width - (glbl_window_width)) / (glbl_window_width / glbl_interval_div));
	if(amount <= 0) {
		amount = floor(width / glbl_window_width);
	}
	if(amount > 15) {
		amount = 15;
	}
	return amount;	
	
}

float** mel_fb(int width, int banks)
{

	float minmel   = 2595*log10(1 + (0/700));
	float maxmel   = 2595*log10(1 + (8000/700));     
	float binwidth = (maxmel - minmel) / (banks + 1);
	float start = minmel; 
	float* mel = calloc(banks + 2, sizeof(float));
	float* h = calloc(banks + 2, sizeof(float)); 
	float* f = calloc(banks + 2, sizeof(float));
	
	float** H = (float**)malloc(width * sizeof(float*));
	//if(H == NULL)
	//	printf("Failed to malloc H in mel_fb\n");
	for(int i = 0; i < width; i++) {
		H[i] = (float*)calloc(banks + 2, sizeof(float));
		//if(H[i] == NULL)
		//	printf("Failed to calloc H[i] in mel_fb\n");
		for(int n = 0; n < banks + 2; n++) {
			H[i][n] = 0;
		}
	}
	
	for(int n = 0; n < banks + 2; n++) {
		mel[n] = start; 
		h[n] = 700*(exp(mel[n]/1125) - 1);
		f[n] = floor(h[n] * (glbl_nfft + 1) / 16000); // 32, 64, 128, 256, 512 (orignally using 128)
		start = start + binwidth;
	}
	
	for(float m = 1; m < banks + 1; m++) {
		for(float k = 1; k < width; k++) {
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
 
	free(mel);
	free(h);
	free(f);
        
	return H;
}

float* mel(float* array, int width, int banks)
{
	float* result = (float*)calloc(banks, sizeof(float));
	float** melfb = mel_fb(width, banks);
	float** applied = calloc(width, sizeof(float*));
	// if(result == NULL || melfb == NULL || applied == NULL)
		// printf("Malloc/Calloc error in mel\n");
	for(int i = 0; i < width; i++) {
		applied[i] = calloc(banks, sizeof(float));
		// if(applied[i] == NULL)
		// printf("Failed to calloc applied[i]\n");
		for(int j = 0; j < banks; j++) {
			applied[i][j] = 0;
		}
	}
	
	for(int i = 0; i < width; i++) {
		for(int j = 0; j < banks; j++) {
			applied[i][j] = array[i] * melfb[i][j];
			if(fpclassify(melfb[i][j]) == FP_INFINITE || fpclassify(melfb[i][j]) == FP_NAN) {
				// printf("MelFb result contained NAN of INF :: %.3f || width := %d || i := %d || j := %d\n", melfb[i][j],  width, i, j);
	 		}
			if(fpclassify(array[i]) == FP_INFINITE || fpclassify(array[i]) == FP_NAN) {
				// printf("Array contained NAN of INF :: %.3f || width := %d\n", melfb[i][j],  width);
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
	for(int i = 0; i < banks; i++) {
		result[i] = result[i] / banks;
	}
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

float* f_realloc(float* input, int length)
{

	float* new = realloc(input, sizeof(float) * length);
	
	if (new != NULL) {
		return new;
	} else {
		fprintf(stderr, "Memory error, trying to allocate %d length which is %d bytes.\n", (int)length, (int)(length*sizeof(float)));
		// free(input);
		exit(-1);
	}
	
}

float* f_paa(float* sequence, int length)
{

	int n_length = floor(length/glbl_paa);
	
	float* temp = malloc(sizeof(float) * length);
	for(int i = 0; i < length; i++) {
		memcpy(&temp[i], &sequence[i], sizeof(float));
	}
	float* new = f_realloc(sequence, n_length);
	
	if (new == NULL) {
		// fprintf(stderr, "Memory error, trying to allocate %d length which is %d bytes.\n", (int)length, (int)(length*sizeof(float)));
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

short* s_realloc(short* input, int length)
{

	short* new = realloc(input, sizeof(short) * length);

	if (new != NULL) {
		return new;
	} else {
		fprintf(stderr, "Memory error, trying to allocate %d length which is %d bytes.\n", (int)length, (int)(length*sizeof(short)));
		exit(-1);
	}
}	

short* paa(short* sequence, int length)
{

	int n_length = floor(length/glbl_paa);
	
	short* temp = malloc(sizeof(short) * length);
	for(int i = 0; i < length; i++) {
		memcpy(&temp[i], &sequence[i], sizeof(short));
	}
	short* new = s_realloc(sequence, n_length);

	if (new == NULL) {
		// fprintf(stderr, "Memory error, trying to allocate %d length which is %d bytes.\n", (int)length, (int)(length*sizeof(short)));
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

float* dct(float* array, int width)
{
	float* result = (float*)calloc(width, sizeof(float));
	if(result == NULL) 
		printf("Failed to malloc 'result' in FFT\n");
	
	for(int i = 0; i < width; i++) {
		float sum = 0;
		for(int j = 0; j < width; j++) {
			sum += array[j] * cos((j + 0.5) * i * (M_PI / width)); 
		}
		result[i] = sum;
	}
	return result;
}

float complex* fft(float complex* chunk, int length)
{
	int n = length; 
	if (n == 1)
		return chunk;

	double complex w = cexp( (2 * M_PI * I) / n);
	float complex* even = calloc(length / 2, sizeof(float complex));
	float complex* odd = calloc(length / 2, sizeof(float complex));
	if(odd == NULL || even == NULL) {
		// printf("Failed to malloc 'odd || even' in FFT\n");
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
		// printf("Failed to malloc 'y' in FFT\n");
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

float* simple_hanning(float* array, int num)
{
	float* results = (float*)calloc(num, sizeof(float));
	float mult = 0;
	for(int i = 0; i < num; i++) {
		mult = 0.5 * (1 - cos(2*M_PI*(i+1)/(num + 1)));
		results[i] = mult * array[i];
	}
	
	return results;
}

float** hanning_chunks_no_accum(float* input, int length, int incr)
{
	int last = floor((length - (incr)) / (incr / glbl_interval_div)) + 1;
	float** chunks = (float**)malloc( (last) * sizeof(float*));
	for(int i = 0; i < (last); i++) {
		chunks[i] = (float*)calloc(incr, sizeof(float));
	}
	int n = 0;
	for(int i = 0; i < (last - 1); i++) {
		float* temp = (float*)malloc(incr * sizeof(float));
		int first = (i) * (incr / glbl_interval_div);
		memcpy(temp, &input[first], incr * sizeof(float));
		float* hanned = simple_hanning(temp, incr);
		memcpy(chunks[n], hanned, sizeof(float) * incr);
		free(hanned);
		n++;
		free(temp);
	}

	return chunks;
}

float** hanning_chunks_no_accum_or_overlap(float* input, int length, int incr)
{
	int last = floor(length / incr);
	float** chunks = (float**)malloc( (last) * sizeof(float*));
	for(int i = 0; i < (last); i++) {
		chunks[i] = (float*)calloc(incr,  sizeof(float));
	}
	int n = 0;
	for(int i = 0; i < (last - 1); i++) {
		float* temp = (float*)malloc(incr * sizeof(float));
		int first = (i) * (incr);
		memcpy(temp, &input[first], incr * sizeof(float));
		float* hanned = simple_hanning(temp, incr);
		memcpy(chunks[n], hanned, sizeof(float) * incr);
		free(hanned);
		n++;
		free(temp);
	}

	return chunks;
}

float* mfcc_quick(float* sequence, int width, int incr, int banks, int paa)
{

	if(paa == 0) {
		sequence = f_paa(sequence, width);
		width = floor(width / glbl_paa);
	}
	int amount = floor((width - (incr)) / (incr / glbl_interval_div)) + 1;

	/* Hanning chunks */
	float** chunks = hanning_chunks_no_accum(sequence, width, incr);
	if(amount <= 0) {
		free(chunks);
		chunks = NULL;
	}
	if(chunks == NULL) {
		chunks = hanning_chunks_no_accum_or_overlap(sequence, width, incr);
		amount = floor(width / incr);
		if(chunks == NULL) {
			return NULL;
		}
	}
	/* Frequency domain */
	float** mags = fft_chunks(chunks, amount, incr);
	for(int i = 0; i < amount; i++) { 
		for (int low = 0, high = (incr / 2) - 1; low < high; low++, high--) {
			float temp = mags[i][low];
			mags[i][low] = mags[i][high];
			mags[i][high] = temp;
		}
	}
	/* Mel filter banks */
	float** applied_mels = malloc(amount * sizeof(float*));
	if(chunks == NULL || mags == NULL || applied_mels == NULL) {
		return NULL;
	}
		
	for(int i = 0; i < amount; i++) {
		applied_mels[i] = (float*)calloc(banks,  sizeof(float));
		if(applied_mels[i] == NULL) {
			return NULL;
		}
		float* temp = mel(mags[i], (incr / 2), banks);
		memcpy(applied_mels[i], temp, sizeof(float) * banks);
		free(temp);
	}
	for(int i = 0; i < amount; i++) {
		for(int j = 0; j < banks; j++) {
			if(applied_mels[i][j] <= 0) {
				applied_mels[i][j] = FLT_EPSILON;
			}
			applied_mels[i][j] = log10(applied_mels[i][j]);	
		}
	}
	for(int i = 0; i < amount; i++) {
		float* temp = dct(applied_mels[i], banks);
		memcpy(applied_mels[i], temp, sizeof(float) * banks);
		free(temp);
	}
	int trunc = floor(banks * glbl_test_trunc);
	float* result = (float*)calloc((amount * trunc), sizeof(float));
	int last = amount;
	int n = 0;
	for(int i = 0; i < amount; i++) {
		for(int j = 0; j < trunc; j++) {
			if(fpclassify(applied_mels[i][j]) == FP_INFINITE || fpclassify(applied_mels[i][j]) == FP_NAN) {
				result[n] = 0; // FLT_EPSILON;
			} else {
				result[n] = applied_mels[i][j];
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
	
	return result;
}
