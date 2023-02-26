#include "knn.h"

int guesscomp(const void * a, const void * b)
{
	long double fa = ((struct Guess*)a)->diff;
	long double fb = ((struct Guess*)b)->diff;
	return (fa > fb) - (fa < fb);
}

int knn_mfccs_size(float* test, int test_length, int k)
{
	int result = 0, n = 0, to_test = 0;
	int mfcc_length = mfcc_size(test_length);
	for(int i = 0; i < num_ph; i++) {
		phones[i] = get_mfcc(phones[i], mfcc_length);
		if(phones[i]->use_count == 0) {
			continue;
		}
		n += phones[i]->use_count;
	}
	if(n == 0) {
		printf("Zero prototypes with size %d, returning sil\n", mfcc_length);
		return (num_ph - 1);
	}
	to_test = n;
	
	struct Guess* gs = (struct Guess*)malloc(sizeof(struct Guess) * to_test);
	int l = 0;
	
	for(int i = 0; i < num_ph; i++) {
		for(int j = 0; j < phones[i]->use_count; j++) {
			gs[l].diff = dtw_frame_result(test, test_length, phones[i], j);
			gs[l].guess = i;
			gs[l].ref_indx = j;
			gs[l].ref = phones[i];
			l++;
		       
		}
	       
	}
		
	qsort(gs, to_test, sizeof(struct Guess), guesscomp);
	int* modes = (int*)calloc(num_ph, sizeof(int));
	
	for(int i = 0; i < num_ph; i++) {
		modes[i] = 0;
	}
	if(k >= to_test) {
		k = floor(to_test / 3);
		if(k <= 0) {
			k = 1;
		}
	}
	for(int i = 0; i < k; i++) {
		modes[gs[i].guess]++;
	}
	int most = 0;
	for(int i = 1; i < num_ph; i++) {
		if(modes[i] > most) {
			most = modes[i];
			result = i;
		}
	}
	free(modes);
	free(gs);
	return result;
}
