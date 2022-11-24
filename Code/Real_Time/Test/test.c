#include "test.h"

struct Phoneme** phones;

char* p_codes[] = {      "b", "d", "k", "p", "t", "g",                                                       // Stops [1 - 6] - OBSTRUENT
	                 "jh", "ch",                                                                         // Affri [7 - 8] - OBSTRUENT
			 "s", "sh", "th", "v", "f", "dh", "z",                                               // Frics [9 - 15] - OBSTRUENT
			 "m", "n", "ng",                                                                     // Nasal [16 - 18] - SONORANT
			 "l", "r", "hh", "w", "y",                                                           // Semiv [19 - 23] - SONORANT
	                 "aa", "ae", "ah", "aw", "er", "ay", "eh", "ey", "ih", "iy", "ow", "oy", "uh", "uw", // Vowel [24 - 37] - SONORANT
			 "sil",                                                                              // Other [38] - OTHER
	           "\0"};

double cubic_interpolate(short y0, short y1, short y2, short y3, double mu)
{
	short a0 = 0, a1 = 0, a2 = 0, a3 = 0, mu2  = 0;

	mu2 = mu*mu;
	a0 = y3 - y2 - y0 + y1;
	a1 = y0 - y1 - a0;
	a2 = y2 - y0;
	a3 = y1;
	return(a0*mu*mu2+a1*mu2+a2*mu+a3);
}

short* resize(short* shorter, size_t s, size_t l)
{
	short* new_sequence = calloc(l, sizeof(short));
	for(unsigned int i = 0; i < s; i++) {
		new_sequence[i] = shorter[i];
	}
	free(shorter);
	size_t i = floor(s / 2);
	while(s < l) {
		for(size_t j = s-1; j >= i ; j--) {
			new_sequence[ j + 1 ] = new_sequence[ j ];
		}
	        		
	        new_sequence[i] = cubic_interpolate((short)new_sequence[i - 2], (short)new_sequence[i - 1], (short)new_sequence[i + 1], (short)new_sequence[i + 2], 0.5);
		s++;
		if(i + 2 >= s - 5){  i = 5; }
		i+=2;	
	}
	return new_sequence;	
}

int test_phoneme_utterance(short* h, int signal_length)
{ 
	int new_length = signal_length;
	if(signal_length % glbl_window_width != 0) {
		while(new_length % glbl_window_width != 0) {
			new_length++;
		}
		h = resize(h, signal_length, new_length);
		signal_length = new_length;
	}
	if((signal_length / glbl_paa) < glbl_window_width) {
		h = resize(h, signal_length, (glbl_window_width * glbl_paa));
		signal_length = (glbl_window_width * glbl_paa);
	}

	float* signal = (float*)malloc(sizeof(float) * signal_length);
	for(int i = 0; i < signal_length; i++) {
		signal[i] = h[i];
	}
	signal = mfcc_quick(signal, signal_length, glbl_window_width, glbl_banks, glbl_paa_op);

	if(signal == NULL) {
		for(int i = 1; i < num_ph; i++) {
			phones[i]->score = DBL_MAX;
		}
		return - 1;
	}
	int new_size = mfcc_size(signal_length);
	if(new_size <= 0) {
		for(int i = 1; i < num_ph; i++) {
			phones[i]->score = DBL_MAX;
		}
		return -1;
	}

	int knn = knn_mfccs_size(signal, signal_length, 7);
	return knn;
}
