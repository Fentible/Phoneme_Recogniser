/**
 * \file knn.c
 * \brief Contains the functions comparing phonemes using KNN
 *
 * Functions for determining the phoneme, phoneme's group and
 * phoneme's voice using MFCCs and raw time domain signals.
 *
 */
   
#include "knn.h"

int knn_mfccs_group(float** test, int test_length, int k, char* ph);
int knn_mfccs_voice(float** test, int test_length, int k, char* ph);

/**
 * \fn guesscomp()
 * \brief A comparator for sorting the results from the distance metric
 */ 
int guesscomp(const void * a, const void * b)
{
	double fa = ((struct Guess*)a)->diff;
	double fb = ((struct Guess*)b)->diff;
	return (fa > fb) - (fa < fb);
}

/** 
 * \fn knn_mfccs()
 * \brief The base KNN function, used to find the final result
 * @test The data which contains the test MFCC, and frame-by-frame zero cross and short time energy arrays
 * @test_length The length of the MFCC as a 1D array
 * @k The amount of guess to be considered in KNN
 * @ph The test phonemes string used only for keeping track of correct and incorrect guesses
 * @prev_ph Used in the \file gram.c grammar functions
 */
int knn_mfccs(float** test, int test_length, int k, char* ph)
{
	time_t started = time(NULL);
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;
	int grp = -1;
	int start = 0, end = 0;
	if(GROUP) {
		grp = knn_mfccs_group(test, test_length, k, ph);
		if      (grp == 0) { start = 1;  end = 7;       }
		else if (grp == 1) { start = 7;  end = 9;       }
		else if (grp == 2) { start = 9;  end = 16;      }
		else if (grp == 3) { start = 16; end = 19;      }
		else if (grp == 4) { start = 19; end = 24;      }
		else if (grp == 5) { start = 24; end = 38;      }
		else               { start = 38; end = 39;      }
	} else {
		start = 1;
		end = num_ph;
	}
	// int n_ph = (end - start);
	int result = 0, n = 0, to_test = 0;
	int mfcc_length = mfcc_size(test_length);
	int n_ph = (end - start);
	for(int i = start; i < end; i++) {
		for(int j = 0; j < phones[i]->size_count; j++) {
				n++;
		}
	}
	if(n < (n_ph)) {
		n = (n_ph);
	}
	
	to_test = n;
	struct Guess* gs = (struct Guess*)malloc(sizeof(struct Guess) * to_test);
	int l = 0;
	
	for(int i = start; i < end; i++) {
		if(to_test == (n_ph)) {
			int sml = INT_MAX, indx = 0;
			for(int j = 0; j < phones[i]->size_count; j++) {
				if(abs(mfcc_length - phones[i]->size[j]) < sml && phones[i]->size[j] != 0) {
					sml = abs(mfcc_length - phones[i]->size[j]);
					indx = j;
				}
			}
			if(phones[i]->size[indx] == 0) {
				printf("Iterating removed all sequences for :: %s :: Ending tests\n", phones[i]->index->name);
				exit(-1);
			}
			gs[l].diff = dtw_frame_result(test[0], test_length, phones[i], indx, glbl_dtw_window);
			gs[l].guess = i;
			gs[l].ref_indx = indx;
			gs[l].ref = phones[i];
			l++;
		} else {
			for(int j = 0; j < phones[i]->size_count; j++) {
				if(phones[i]->size[j] != 0) {
					gs[l].diff = dtw_frame_result(test[0], test_length, phones[i], j, glbl_dtw_window);
					gs[l].guess = i;
					gs[l].ref_indx = j;
					gs[l].ref = phones[i];
					l++;
				}
			}
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
		if(strcmp(gs[i].ref->index->name, ph) != 0) {
			gs[i].ref->error[gs[i].ref_indx]++;
		} else {
			gs[i].ref->correct[gs[i].ref_indx]++;
		}
	}
	int most = 0, final = 0;
	for(int i = 1; i < num_ph; i++) {
		if(modes[i] > most) {
			most = modes[i];
			final = i;
		}
	}
	free(gs);
	free(modes);
	result = final;
	time_t ended = time(NULL);
	total_test_time += (ended - started);
	total_dtw_tests++;
	return result;
}

/** 
 * \fn knn_mfccs_size()
 * \brief The base KNN function, used to find the final result
 * @total_test The data which contains the test MFCC, and frame-by-frame zero cross and short time energy arrays
 * @test_length The length of the MFCC as a 1D array
 * @k The amount of guess to be considered in KNN
 * @ph The test phonemes string used only for keeping track of correct and incorrect guesses
 * @prev_ph Used in the \file gram.c grammar functions
 * @return The guessed phoneme's index 
 *
 * Unlike \fn knn_mfccs() this version only compares those of the same size
 */
int knn_mfccs_size(float** total_test, int test_length, int k, char* ph)
{
	time_t started = time(NULL);
	float* test = total_test[0];
	int j = 1;
	//int trunc = floor(glbl_banks * glbl_test_trunc);
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;
	
	int grp = -1;
	int start = 0, end = 0;
	if(GROUP) {
		// Quick and easy ~50%
		/* ZC = 1; */
		/* STE = 0; */
		/* grp = knn_mfccs_group(total_test, test_length, k, ph); */
		/* STE = 0; */
		/* ZC = 0; */
		/* for(int i = 1; i < num_ph; i++) { */
		/* 	if(strcmp(ph, phones[i]->index->name) == 0) { */
		/* 		grp = phones[i]->index->group_i; */
		/* 	} */
		/* } */

		// Combined method
		// int* zc_group = NULL;
		// int* ste_group = NULL;
		// int* kurt_group = NULL;
		// int* entr_group = NULL;
		
		STE = 1;
		// ste_group = knn_mfccs_group_comb(total_test, test_length, k, ph);
		STE = 0;
		
		ZC = 1;
		// zc_group = knn_mfccs_group_comb(total_test, test_length, k, ph);
		ZC = 0;

		KURT = 1;
		// kurt_group = knn_mfccs_group_comb(total_test, test_length, k, ph);
		KURT = 0;

		ENTR = 1;
		// entr_group = knn_mfccs_group_comb(total_test, test_length, k, ph);
		ENTR = 0;
		
		/* int lrgst = 0; */
		/* for(int i = 0; i < 7; i++) { */
		/* 	if((zc_group[i]) > lrgst) { */
		/* 		lrgst = zc_group[i]; */
		/* 		grp = i; */
		/* 	} */
		/* } */

		grp = knn_mfccs_group(total_test, test_length, k, ph);

		// free(ste_group);
		// free(zc_group);
		// free(kurt_group);
		// free(entr_group);

		if(grp == -1) {
			printf("Whoops... group is -1 \n");
			exit(-1);
		}
		if      (grp == 0) { start = 1;  end = 7;       }
		else if (grp == 1) { start = 7;  end = 9;       }
		else if (grp == 2) { start = 9;  end = 16;      }
		else if (grp == 3) { start = 16; end = 19;      }
		else if (grp == 4) { start = 19; end = 24;      }
		else if (grp == 5) { start = 24; end = 38;      }
		else               { return (num_ph - 1);       } // return num_ph - 1?
	} else if (VOICED && !(GROUP)) {
		ZC = 1;
		grp = knn_mfccs_voice(total_test, test_length, k, ph);
		ZC = 0;
		if(grp == -1) {
			printf("Whoops... \n");
			exit(-1);
		}
		if      (grp == 0) { start = 1;  end = 16; }
		else if (grp == 1) { start = 16; end = 38; }
		else               { return (num_ph - 1); } // return 6?
	} else {
		start = 1;
		end = num_ph;
	}
	// int n_ph = (end - start);
	int result = 0, n = 0, to_test = 0;
	int mfcc_length = mfcc_size(test_length);
	int n_ph = (end - start);
	for(int i = start; i < end; i++) {
		for(int j = 0; j < phones[i]->size_count; j++) {
			if(phones[i]->size[j] != 0 && mfcc_length == phones[i]->size[j]) {
				n++;
			}
		}
	}
	int NONE = 0;
	if(n == 0) {
		n = (n_ph);
		NONE = 1;
	}
	
	to_test = n;
	struct Guess* gs = (struct Guess*)malloc(sizeof(struct Guess) * to_test);
	int l = 0;
	
	for(int i = start; i < end; i++) {
		if(NONE == 1) {
			int sml = INT_MAX, indx = 0;
			for(int j = 0; j < phones[i]->size_count; j++) {
				if(abs(mfcc_length - phones[i]->size[j]) < sml && phones[i]->size[j] != 0) {
					sml = abs(mfcc_length - phones[i]->size[j]);
					indx = j;
				}
			}
			if(phones[i]->size[indx] == 0) {
				printf("Iterating removed all sequences for :: %s :: Ending tests\n", phones[i]->index->name);
				exit(-1);
			}
			gs[l].diff = dtw_frame_result(test, test_length, phones[i], indx, glbl_dtw_window);
			gs[l].guess = i;
			gs[l].ref_indx = indx;
			gs[l].ref = phones[i];
			l++;
		} else {
			for(int j = 0; j < phones[i]->size_count; j++) {
				if(phones[i]->size[j] != 0 && mfcc_length == phones[i]->size[j]) {
					gs[l].diff = dtw_frame_result(test, test_length, phones[i], j, glbl_dtw_window);
					gs[l].guess = i;
					gs[l].ref_indx = j;
					gs[l].ref = phones[i];
					l++;
				}
			}
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
		if(ph == NULL)
			continue;
		if(strcmp(gs[i].ref->index->name, ph) != 0) {
			gs[i].ref->error[gs[i].ref_indx]++;
		} else {
			gs[i].ref->correct[gs[i].ref_indx]++;
		}
	}
	int most = 0, final = 0;
	for(int i = 1; i < num_ph; i++) {
		if(modes[i] > most) {
			most = modes[i];
			final = i;
		}
	}
	free(gs);
	free(modes);

	result = final;
	
	
	time_t ended = time(NULL);
	total_test_time += (ended - started);
	total_dtw_tests++;
	return result;
}

/** 
 * \fn knn_mfccs_group()
 * \brief The KNN function used to determine the phoneme group
 * @test The data which contains the test MFCC, and frame-by-frame zero cross and short time energy arrays
 * @test_length The length of the MFCC as a 1D array
 * @k The amount of guess to be considered in KNN
 * @ph The test phonemes string used only for keeping track of correct and incorrect guesses
 * @prev_ph Used in the \file gram.c grammar functions
 * @return The guessed group to be used in a base KNN function such as \fn knn_mfcc_size()
 *
 * This function decides which group a phoneme is in i.e. stop, nasal, etc. and the result is to be handled in another function
 */
int knn_mfccs_group(float** test, int test_length, int k, char* ph)
{
	k = glbl_group_k;
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;

	int grp = -1;
	int start = 0, end = 0;
	if(VOICED) {
		// int* zc_group = NULL;
		// int* ste_group = NULL;
		/* STE = 1; */
		/* ste_group = knn_mfccs_voice_comb(test, test_length, k, ph); */
		/* STE = 0; */
		
		/* ZC = 1; */
		/* zc_group = knn_mfccs_voice_comb(test, test_length, k, ph); */
		/* ZC = 0; */
		/* int lrgst = 0; */
		/* for(int i = 0; i < 3; i++) { */
		/* 	if((ste_group[i] + zc_group[i]) > lrgst) { */
		/* 		lrgst = ste_group[i] + zc_group[i]; */
		/* 		grp = i; */
		/* 	} */
		/* } */
		/* free(ste_group); */
		/* free(zc_group); */
		ZC = 1;
		grp = knn_mfccs_voice(test, test_length, k, ph);
		ZC = 0;
		if(grp == -1) {
			printf("Whoops... \n");
			exit(-1);
		}
		if      (grp == 0) { start = 1;  end = 16; }
		else if (grp == 1) { start = 16; end = 38; }
		else               { return 6; } // return 6?
	} else {
		start = 1;
		end = num_ph;
	}

	int result = 0, n = 0, to_test = 0;
	int mfcc_length = mfcc_size(test_length);
	int n_ph = (end - start);
	for(int i = start; i < end; i++) {
		for(int j = 0; j < phones[i]->size_count; j++) {
			if(mfcc_length == phones[i]->size[j]) {
				n++;
				// break;
			}
		}

	}
	if(n < (n_ph)) {
		n = (n_ph);
	}
	to_test = n;
	struct Guess* gs = (struct Guess*)malloc(sizeof(struct Guess) * to_test);
	int l = 0;
	
	for(int i = start; i < end; i++) {
		if(to_test == n_ph) {
			int sml = INT_MAX, indx = 0;
			for(int j = 0; j < phones[i]->size_count; j++) {
				if(abs(mfcc_length - phones[i]->size[j]) < sml && phones[i]->size[j] != 0) {
					sml = abs(mfcc_length - phones[i]->size[j]);
					indx = j;
				}
			}
			if(STE || ZC || DELTA || DELTA_DELTA) {
				gs[l].diff = dtw_frame_result_group(test, test_length, phones[i], indx, glbl_dtw_window);
			} else {
				gs[l].diff = dtw_frame_result(test[0], test_length, phones[i], indx, glbl_dtw_window);
			}
			gs[l].guess = i;
			l++;
		} else {
			for(int j = 0; j < phones[i]->size_count; j++) {
				if(mfcc_length == phones[i]->size[j]) {
					if(STE || ZC || DELTA || DELTA_DELTA) {
						gs[l].diff = dtw_frame_result_group(test, test_length, phones[i], j, glbl_dtw_window);
					} else {
						gs[l].diff = dtw_frame_result(test[0], test_length, phones[i], j, glbl_dtw_window);
					}
					gs[l].guess = i;
					l++;
					// break;
				}
			}
		}

	}
	// {"STOP", "AFRI", "FRIC", "NASL", "SEMV", "VOWL", "OTHR"};

	qsort(gs, to_test, sizeof(struct Guess), guesscomp);
	int* group = (int*)calloc(7, sizeof(int));
	
	for(int i = 0; i < k; i++) {
		// modes[gs[i].guess]++;
		if      (gs[i].guess >=  1 && gs[i].guess <=  6) { group[0]++;  }
		else if (gs[i].guess ==  7 || gs[i].guess ==  8) { group[1]++;  }
		else if (gs[i].guess >=  9 && gs[i].guess <= 15) { group[2]++;  }
		else if (gs[i].guess >= 16 && gs[i].guess <= 18) { group[3]++;  }
		else if (gs[i].guess >= 19 && gs[i].guess <= 23) { group[4]++;  }
		else if (gs[i].guess >= 24 && gs[i].guess <= 37) { group[5]++;  }
		else                                             { group[6]++;  }
	}
	int most = 0, final = 0;
	for(int i = 0; i < 7; i++) {
		if(group[i] > most) {
			most = group[i];
			final = i;
		}
	}
	if(!BOUNDS) { 
		int indx = 0;
		for(int i = 1; i < num_ph; i++) {
			if(strcmp(ph, phones[i]->index->name) == 0) {
				indx = i;
			}
		}

		group_matrix[phones[indx]->index->group_i][final]++;
	}
	
	free(gs);
	free(group);
	result = final;
	return result;
}

/** 
 * \fn knn_mfccs_voice()
 * \brief The KNN function used to determine the phoneme voice type
 * @test The data which contains the test MFCC, and frame-by-frame zero cross and short time energy arrays
 * @test_length The length of the MFCC as a 1D array
 * @k The amount of guess to be considered in KNN
 * @ph The test phonemes string used only for keeping track of correct and incorrect guesses
 * @prev_ph Used in the \file gram.c grammar functions
 * @return The guessed voice type to be used in a group KNN function such as \fn knn_mfcc_group()
 *
 * This function decides whether a test sequence is voiced, unvoiced or silence and the result is to be handled in another function
 */
int knn_mfccs_voice(float** test, int test_length, int k, char* ph)
{
	k = glbl_voice_k;
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;
	int result = 0, n = 0, to_test = 0;
	int mfcc_length = mfcc_size(test_length);
	for(int i = 1; i < num_ph; i++) {
		for(int j = 0; j < phones[i]->size_count; j++) {
			if(mfcc_length == phones[i]->size[j]) {
				n++;
				// break;
			}
		}

	}
	if(n < (num_ph - 1)) {
		n = (num_ph - 1);
	}
	to_test = n;
	struct Guess* gs = (struct Guess*)malloc(sizeof(struct Guess) * to_test);
	int l = 0;
	
	for(int i = 1; i < num_ph; i++) {
		if(to_test == (num_ph - 1)) {
			int sml = INT_MAX, indx = 0;
			for(int j = 0; j < phones[i]->size_count; j++) {
				if(abs(mfcc_length - phones[i]->size[j]) < sml && phones[i]->size[j] != 0) {
					sml = abs(mfcc_length - phones[i]->size[j]);
					indx = j;
				}
			}
			if(STE || ZC || DELTA || DELTA_DELTA) {
				gs[l].diff = dtw_frame_result_group(test, test_length, phones[i], indx, glbl_dtw_window);
			} else {
				gs[l].diff = dtw_frame_result(test[0], test_length, phones[i], indx, glbl_dtw_window);
			}
			gs[l].guess = i;
			l++;
		} else {
			for(int j = 0; j < phones[i]->size_count; j++) {
				if(mfcc_length == phones[i]->size[j]) {
					if(STE || ZC || DELTA || DELTA_DELTA) {
						gs[l].diff = dtw_frame_result_group(test, test_length, phones[i], j, glbl_dtw_window);
					} else {
						gs[l].diff = dtw_frame_result(test[0], test_length, phones[i], j, glbl_dtw_window);
					}
					gs[l].guess = i;
					l++;
					// break;
				}
			}
		}

	}
	// {"STOP", "AFRI", "FRIC", "NASL", "SEMV", "VOWL", "OTHR"};

	qsort(gs, to_test, sizeof(struct Guess), guesscomp);
	int* group = (int*)calloc(3, sizeof(int));
	
	for(int i = 0; i < k; i++) {
		// modes[gs[i].guess]++;
		if      (gs[i].guess >=  1 && gs[i].guess <=  6) { group[0]++;  }
		else if (gs[i].guess ==  7 || gs[i].guess ==  8) { group[0]++;  }
		else if (gs[i].guess >=  9 && gs[i].guess <= 15) { group[0]++;  }
		else if (gs[i].guess >= 16 && gs[i].guess <= 18) { group[1]++;  }
		else if (gs[i].guess >= 19 && gs[i].guess <= 23) { group[1]++;  }
		else if (gs[i].guess >= 24 && gs[i].guess <= 37) { group[1]++;  }
		else                                             { group[2]++;  }
	}
	int most = 0, final = 0;
	for(int i = 0; i < 3; i++) {
		if(group[i] > most) {
			most = group[i];
			final = i;
		}
	}
	if(!BOUNDS) { 
		int indx = 0;
		for(int i = 1; i < num_ph; i++) {
			if(strcmp(ph, phones[i]->index->name) == 0) {
				indx = i;
			}
		}
	
		voice_matrix[phones[indx]->index->voice][final]++;
	}
	
	free(gs);
	free(group);
	result = final;
	return result;
}

int knn_mfccs_voice_time(float* test, int test_length, int k, char* ph)
{
	k = glbl_group_k;
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;
	int result = 0, n = 0, to_test = 0;
	int mfcc_length = test_length;
	for(int i = 1; i < num_ph; i++) {
		for(int j = 0; j < phones[i]->raw_count; j++) {
			// if(mfcc_length == phones[i]->raw_sizes[j]) {
				n++;
				// break;
				// }
		}

	}
	if(n < (num_ph - 1)) {
		n = (num_ph - 1);
	}
	to_test = n;
	struct Guess* gs = (struct Guess*)malloc(sizeof(struct Guess) * to_test);
	int l = 0;
	
	for(int i = 1; i < num_ph; i++) {
		if(to_test == (num_ph - 1)) {
			int sml = INT_MAX, indx = 0;
			for(int j = 0; j < phones[i]->raw_count; j++) {
				if(abs(mfcc_length - phones[i]->raw_sizes[j]) < sml) {
					sml = abs(mfcc_length - phones[i]->raw_sizes[j]);
					indx = j;
				}
			}
			gs[l].diff = dtw_frame_result_group_time(test, test_length, phones[i], indx, glbl_dtw_window);
			gs[l].guess = i;
			l++;
		} else {
			for(int j = 0; j < phones[i]->raw_count; j++) {
				// if(mfcc_length == phones[i]->raw_sizes[j]) {
					gs[l].diff = dtw_frame_result_group_time(test, test_length, phones[i], j, glbl_dtw_window);
					gs[l].guess = i;
					l++;
					// break;
					// }
			}
		}

	}
	// {"STOP", "AFRI", "FRIC", "NASL", "SEMV", "VOWL", "OTHR"};

	qsort(gs, to_test, sizeof(struct Guess), guesscomp);
	int* group = (int*)calloc(3, sizeof(int));
	
	for(int i = 0; i < k; i++) {
		// modes[gs[i].guess]++;
		if      (gs[i].guess >=  1 && gs[i].guess <=  6) { group[0]++;  }
		else if (gs[i].guess ==  7 || gs[i].guess ==  8) { group[0]++;  }
		else if (gs[i].guess >=  9 && gs[i].guess <= 15) { group[0]++;  }
		else if (gs[i].guess >= 16 && gs[i].guess <= 18) { group[1]++;  }
		else if (gs[i].guess >= 19 && gs[i].guess <= 23) { group[1]++;  }
		else if (gs[i].guess >= 24 && gs[i].guess <= 37) { group[1]++;  }
		else                                             { group[2]++;  }
	}
	int most = 0, final = 0;
	for(int i = 0; i < 3; i++) {
		if(group[i] > most) {
			most = group[i];
			final = i;
		}
	}
	if(!BOUNDS && ph != NULL) { 
		int indx = 0;
		for(int i = 1; i < num_ph; i++) {
			if(strcmp(ph, phones[i]->index->name) == 0) {
				indx = i;
			}
		}
	
		voice_matrix[phones[indx]->index->voice][final]++;
	}
	
	free(gs);
	free(group);
	result = final;
	return result;
}

/** 
 * \fn knn_mfccs_group_time()
 * \brief The KNN function used to determine the phoneme group
 * @test The data which contains the test time domain signal
 * @test_length The length of the MFCC as a 1D array
 * @k The amount of guess to be considered in KNN
 * @ph The test phonemes string used only for keeping track of correct and incorrect guesses
 * @prev_ph Used in the \file gram.c grammar functions
 * @return The guessed voice type to be used in a group KNN function such as \fn knn_mfcc_size()
 *
 * This function decides whether a test sequence is voiced, unvoiced or silence and the result is to be handled in another function
 * Unlike \fn knn_mfccs_group() this function uses a time domain signal only to find out the group.
 */
int knn_mfccs_group_time(float* test, int test_length, int k, char* ph)
{
	k = glbl_group_k;
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;
	int result = 0, n = 0, to_test = 0;
	int mfcc_length = test_length;
	for(int i = 1; i < num_ph; i++) {
		for(int j = 0; j < phones[i]->raw_count; j++) {
			if(mfcc_length == phones[i]->raw_sizes[j]) {
				n++;
				// break;
			}
		}

	}
	if(n < (num_ph - 1)) {
		n = (num_ph - 1);
	}
	to_test = n;
	struct Guess* gs = (struct Guess*)malloc(sizeof(struct Guess) * to_test);
	int l = 0;
	
	for(int i = 1; i < num_ph; i++) {
		if(to_test == (num_ph - 1)) {
			int sml = INT_MAX, indx = 0;
			for(int j = 0; j < phones[i]->raw_count; j++) {
				if(abs(mfcc_length - phones[i]->raw_sizes[j]) < sml) {
					sml = abs(mfcc_length - phones[i]->raw_sizes[j]);
					indx = j;
				}
			}
			gs[l].diff = dtw_frame_result_group_time(test, test_length, phones[i], indx, glbl_dtw_window);
			gs[l].guess = i;
			l++;
		} else {
			for(int j = 0; j < phones[i]->raw_count; j++) {
				if(mfcc_length == phones[i]->raw_sizes[j]) {
					gs[l].diff = dtw_frame_result_group_time(test, test_length, phones[i], j, glbl_dtw_window);
					gs[l].guess = i;
					l++;
					// break;
				}
			}
		}

	}
	// {"STOP", "AFRI", "FRIC", "NASL", "SEMV", "VOWL", "OTHR"};
	qsort(gs, to_test, sizeof(struct Guess), guesscomp);
	int* group = (int*)calloc(7, sizeof(int));
	
	for(int i = 0; i < k; i++) {
		// modes[gs[i].guess]++;
		if      (gs[i].guess >=  1 && gs[i].guess <=  6) { group[0]++;  }
		else if (gs[i].guess ==  7 || gs[i].guess ==  8) { group[1]++;  }
		else if (gs[i].guess >=  9 && gs[i].guess <= 15) { group[2]++;  }
		else if (gs[i].guess >= 16 && gs[i].guess <= 18) { group[3]++;  }
		else if (gs[i].guess >= 19 && gs[i].guess <= 23) { group[4]++;  }
		else if (gs[i].guess >= 24 && gs[i].guess <= 37) { group[5]++;  }
		else                                             { group[6]++;  }
	}
	int most = 0, final = 0;
	for(int i = 0; i < 7; i++) {
		if(group[i] > most) {
			most = group[i];
			final = i;
		}
	}
	if(!BOUNDS && ph != NULL) { 
		int indx = 0;
		for(int i = 1; i < num_ph; i++) {
			if(strcmp(ph, phones[i]->index->name) == 0) {
				indx = i;
			}
		}

		group_matrix[phones[indx]->index->group_i][final]++;
	}
	
	free(gs);
	free(group);
	result = final;
	return result;
}

/** 
 * \fn knn_mfccs_voice_time()
 * \brief The KNN function used to determine the phoneme group
 * @test The data which contains the test time domain signal
 * @test_length The length of the MFCC as a 1D array
 * @k The amount of guess to be considered in KNN
 * @ph The test phonemes string used only for keeping track of correct and incorrect guesses
 * @prev_ph Used in the \file gram.c grammar functions
 * @return The guessed voice type to be used in a group KNN function such as \fn knn_mfcc_group()
 *
 * This function decides whether a test sequence is voiced, unvoiced or silence and the result is to be handled in another function.
 * Unlike \fn knn_mfccs_voice() this function uses a time domain signal only to find out the voice type.
 */
int k_means(float* test, int test_length, int k)
{
	time_t started = time(NULL);

	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;
	int result = 0, n = 0, to_test = 0;
	// int trunc = floor(glbl_banks * glbl_test_trunc);
	// int mfcc_len = mfcc_size(test_length);
	for(int i = 1; i < num_ph; i++) {
		n++;
	}
	to_test = n;

	struct Guess* gs = (struct Guess*)malloc(sizeof(struct Guess) * to_test);
	for(int i = 0; i < to_test; i++) {
		gs[i].diff = INT_MAX;
	}
	int l = 0;
	for(int i = 1; i < num_ph; i++) {
		gs[l].diff = dtw_clust_result(test, test_length, phones[i], 0, 0, glbl_dtw_window);
		gs[l].guess = i;
		l++;
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
	int most = 0, final = 0;
	for(int i = 1; i < num_ph; i++) {
		if(modes[i] > most) {
			most = modes[i];
			final = i;
		}
	}

	free(gs);
	free(modes);
	time_t ended = time(NULL);
	total_test_time += (ended - started);
	total_dtw_tests++;
	result = final;
	return result;
}
