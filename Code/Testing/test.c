/**
 * @file   test.c
 * @author T. Buckingham
 * @date   Wed May  4 20:20:52 2022
 * 
 * @brief  Contains the  functions for classifying test data with the trained model
 * 
 * 
 */

#include "test.h"

struct guess {
	long start;
	long end;
	char* name;
};  /* Structure to store the guess made during testing with boundaries */

int** matrix;        /* The final confusion matrix for phonemes */
float* per_correct;  /* Percentage correct for each phoneme */

int group_matrix[7][7];
int voice_matrix[3][3];
int sil_v_matrix[3][3];

int sil_corr;

int* shorted;
int* removed;

// Correct, incorrect, correct group classifications, amount of tests completed and the amount of classifications that result in more than one with a score of 0
int correct = 0, fails = 0, group = 0, done = 0, zero_fails = 0;

long double min_edit = 0;         /* The minimum edit distance for the whole boundary test cycle */
long double ins = 0;              /* The total amount of insertions for the whole boundary test cycle */
long double delts = 0;            /* The total amount of deletions for the whole boundary test cycle */
long double subs = 0;             /* The total amount of substitutions for the whole boundary test cycle */
int tested_files = 0;             /* The total amount of files tested for the whole test cycle */
long double corr = 0;             /* The total amount of correct classifications for the whole boundary test cycle */
long double WER = 0;              /* The overall WER for the whole boundary test cycle */
long double refs = 0;             /* The total amount of phonemes in the rerference files */
long double avg_wer = 0.0;        /* The average WER for the whole boundary test cycle */
long double low_wer = 0.0;        /* The lowest WER for the whole boundary test cycle */
long double high_wer = 0;         /* The highest WER for the whole boundary test cycle */


long double worst, total_lengths = 0, total_dists = 0, total_edits = 0, worst_case = 0;

int iteration_limit = 0;          /* The maximum amount of test files per iteration */
int current_chunk = 0;            /* The current index of chunks when using SPLIT_DATA */
 
int group_matrix[7][7] = {0};     /* The final confusion matrix for phoneme groups */
int voice_matrix[3][3] = {0};     /* The final confusion matrix for phoneme voice types */

char* best_file;                  /* The filename of the reference file with the lowest WER */
char* worst_file;                 /* The filename of the reference file with the highest WER */

static char* male_dir = "../Testing/MALE/";
static char* female_dir = "../Testing/FEMALE/";
static char* test_dir = "../Testing/TEST/";
static char* spkr1_dir = "../Testing/SPKR1/";
static char* spkr1_nosil_dir = "../Testing/SPKR1_NOSIL/";
static char* spkr1_cafe_dir = "../Testing/SPKR1_CAFE/";
static char* spkr1_car_dir = "../Testing/SPKR1_CAR/";
static char* spkr1_white_dir = "../Testing/SPKR1_WHITE/";
static char* spkr1_strt_dir = "../Testing/SPKR1_STRT/";

static char* res_male_dir = "../Testing/MALE/res/";
static char* res_female_dir = "../Testing/FEMALE/res/";
static char* res_test_dir = "../Testing/TEST/res/";
static char* res_spkr1_dir = "../Testing/SPKR1/res/";
static char* res_spkr1_nosil_dir = "../Testing/SPKR1_NOSIL/res/";
static char* res_spkr1_cafe_dir = "../Testing/SPKR1_CAFE/res/";
static char* res_spkr1_car_dir = "../Testing/SPKR1_CAR/res/";
static char* res_spkr1_white_dir = "../Testing/SPKR1_WHITE/res/";
static char* res_spkr1_strt_dir = "../Testing/SPKR1_STRT/res/";

struct codes {
	char** codes;
	int amount;
}; /* Stores the codes found in a reference files or result (.res) files */

void export_results_pca(char* ph_code);
int test_phoneme_utterance(short* h, int signal_length);
void minimum_edit_distance(char* res_filename, char* phn_filename);
void merge_silences(char* filename);
struct codes* get_codes(char* filename);
void reset(void);

/** 
 * @brief Resets all values used for determining accuracy and outputting results.
 *
 * This function's primary use is to reset the values when SPLIT_DATA is in use as the accuracy values and confusions matrices should be per test loop 
 */
void reset(void)
{
	per_correct = (float*)calloc(num_ph, sizeof(float));
	matrix = (int**)malloc(sizeof(int*) * num_ph);
	for(int i = 1; i < num_ph; i++) {
		matrix[i] = (int*)calloc(num_ph, sizeof(int));
	}
	correct = 0, fails = 0, group = 0, done = 0, zero_fails = 0;
	min_edit = 0, total_lengths = 0, total_dists = 0;
	ins = 0, delts = 0, subs = 0;
	total_edits = 0, tested_files = 0, worst_case = 0;
	corr = 0, worst = 0, WER = 0, refs = 0, avg_wer = 0;
	for(int i = 0; i < 7; i++) {
		for(int j = 0; j < 7; j++) {
			group_matrix[i][j] = 0;
		}
	}
	for(int i = 0; i < 3; i++) {
		for(int j = 0; j < 3; j++) {
			voice_matrix[i][j] = 0;
			sil_v_matrix[i][j] = 0;
		}
	}

	DIR* theFolder;
	char* path;
	if(FEMALE) {
		theFolder = opendir("../Testing/FEMALE/res/");
		path = res_female_dir;
	} else if(MALE) {
		theFolder = opendir("../Testing/MALE/res/");
		path = res_male_dir;
	} else if(SPKR1 && !(CAFE || STRT || WHITE || CAR)) {
		theFolder = opendir("../Testing/SPKR1/res/");
		path = res_male_dir;
	} else if(SPKR1_NOSIL) {
		theFolder = opendir("../Testing/SPKR1_NOSIL/res/");
		path = res_male_dir;
	} else if(SPKR1 && CAFE) {
		theFolder = opendir("../Testing/SPKR1_CAFE/res/");
		path = res_male_dir;
	} else if(SPKR1 && STRT) {
		theFolder = opendir("../Testing/SPKR1_STRT/res/");
		path = res_male_dir;
	} else if(SPKR1 && WHITE) {
		theFolder = opendir("../Testing/SPKR1_WHITE/res/");
		path = res_male_dir;
	} else if(SPKR1 && CAR) {
		theFolder = opendir("../Testing/SPKR1_CAR/res/");
		path = res_male_dir;
	} else {
		theFolder = opendir("../Testing/TEST/res/");
		path = res_test_dir;
	}
	struct dirent *next_file;
	char filepath[512];

	while ( (next_file = readdir(theFolder)) != NULL )
	{
		
		// build the path for each file in the folder
		sprintf(filepath, "%s/%s", path, next_file->d_name);
		remove(filepath);
	}
	closedir(theFolder);
	
	return;
}

/** 
  * @brief Extracts the labels (start, end and code) from a given file
  * 
  * @param filename The name of the filename to extract the labels from
  * 
  * @return A codes structure containing a list of phoneme codes in the same order as written in the file.
  */
struct codes* get_codes(char* filename)
{
	struct codes* res = (struct codes*)malloc(sizeof(struct codes));
	char* end_ptr;
	FILE* fp = fopen(filename, "r");
	if(!fp) { printf("Error opening '%s' file when getting codes\n", filename); exit(-1);}
	char line[256];
	char *p;
	int amount = 0;
	while(fgets(line, sizeof(line), fp)) {
		amount++;
		// printf("%d : %s\n", amount, line);
	}
	fseek(fp, 0, SEEK_SET);
	// amount--; // removing new line
	if(amount == 0) {
		res->amount = 0;
		res->codes = NULL;
		return res;
	}
	res->codes = (char**)malloc(amount * sizeof(char*));
	int n = 0;
	for(int i = 0; i < amount; i++) {
		fgets(line, sizeof(line), fp);
		p = strtok(line," ");
		strtol(p, &end_ptr, 10);
	      	p = strtok(NULL, " ");
		strtol(p, &end_ptr, 10);
		p = strtok(NULL, " ");
		while(p[strlen(p) - 1] == '\n' || p[strlen(p) - 1] ==  ' ') {
			p[strlen(p) - 1] = '\0';
		}
		
		res->codes[n] = (char*)malloc((strlen(p) + 1) * sizeof(char));
		strcpy(res->codes[n], p);
		n++;
	}
	res->amount = amount;
	fclose(fp);
	
	return res;
}

/** 
 * @brief Determines the minimum edit distance between a reference and result file
 * 
 * @param res_filename The result file produced by the testing process
 * @param phn_filename The reference file given by the user in the test data
 */
void minimum_edit_distance(char* res_filename, char* phn_filename)
{

	struct codes* res_codes = get_codes(res_filename);
	struct codes* phn_codes = get_codes(phn_filename);
	// printf("%s : %s\n", res_filename, phn_filename);
	if(res_codes->amount <= 0 || phn_codes->amount <= 0) {
		if(res_codes->codes != NULL) {
			free(res_codes->codes);
		}
		if(phn_codes->codes != NULL) {
			free(phn_codes->codes);
		}
		free(res_codes);
		free(phn_codes);
		delts += phn_codes->amount;
		printf("None found for %s, exiting...\n", res_filename);
		return;
		
	}
	int** matrix = (int**)malloc(sizeof(int*) * (res_codes->amount + 1));
	for(int i = 0; i < res_codes->amount + 1; i++) {
		matrix[i] = (int*)calloc(phn_codes->amount + 1, sizeof(int));
	}
	matrix[0][0] = 0;

	int cost = 0;
	for(int i = 1; i <= res_codes->amount; i++) {
		for(int j = 1; j <= phn_codes->amount; j++) {
			if(strcmp(res_codes->codes[i-1], phn_codes->codes[j-1]) == 0) {
				cost = 0;
			} else {
				cost = 1;
			}
			int temp_min = min((matrix[i - 1][j] + 1), (matrix[i][j - 1] + 1));
			matrix[i][j] = min((matrix[i-1][j-1] + cost), temp_min);
		}
	}

	int row = res_codes->amount;
	int col = phn_codes->amount;
	long double l_subs = subs, l_delts = delts, l_ins = ins;
	while(col != 0 && row != 0) {
		int prev_cost = matrix[row][col];
		int subst = matrix[row-1][col-1];
		int delet = matrix[row-1][col];
		int inser = matrix[row][col-1];
		int temp_min = min(inser, delet);
		int min_cost = min(subst, temp_min);
		if(min_cost == prev_cost) {
			corr++;
			row = row - 1;
			col = col - 1;
		} else if(row != 0 && col != 0 && min_cost == subst) {
			subs++;
			row = row - 1;
			col = col - 1;
		} else if(row != 0 && min_cost == delet) {
			delts++;
			row = row - 1;
			col = col;
		} else if(col != 0 && min_cost == inser) {
			ins++;
			row = row;
			col = col - 1;
		} else {
			printf("Backtrace failed\n");
		}
		total_dists++;
									  
	}
	int dis = matrix[res_codes->amount][phn_codes->amount];
	int larger = max(res_codes->amount, phn_codes->amount);
	total_lengths += abs(res_codes->amount - phn_codes->amount);
	worst += larger;
	min_edit += ((((float)larger -  (float)dis) / (double)larger)) * 100;
	refs += phn_codes->amount;
	float curr_wer = (((ins-l_ins)+(subs-l_subs)+(delts-l_delts))/phn_codes->amount); 
	avg_wer += curr_wer;
	if(curr_wer > high_wer) {
		high_wer = curr_wer;
	        worst_file = (char*)realloc(worst_file, sizeof(char) * (strlen(res_filename) + 1));
		strcpy(worst_file, res_filename);
	}
	if((corr / phn_codes->amount) > low_wer) {
		low_wer = (corr / phn_codes->amount);
		best_file = (char*)realloc(best_file, sizeof(char) * (strlen(res_filename) + 1));
		strcpy(best_file, res_filename);
	}

	for(int i = 0; i < res_codes->amount + 1; i++) {
		free(matrix[i]);
	}
	free(matrix);
	for(int i = 0; i < res_codes->amount; i++) {
		free(res_codes->codes[i]);
	}
	for(int i = 0; i < phn_codes->amount; i++) {
		free(phn_codes->codes[i]);
	}
	free(res_codes->codes);
	free(phn_codes->codes);
	free(res_codes);
	free(phn_codes);
	
	return;
}

/** 
 * @brief Tests files without knowing the boundaries, will use the boundary detection to find them.
 * 
 * @param h The audio signal to test
 * @param signal_length The length of @param h
 * 
 * @return The classified phoneme's index
 */
int test_phoneme_utterance(short* h, int signal_length)
{
	// int test_length = signal_length;
	int sil = is_sil_ste(h, signal_length);
	// int sil_zc = is_sil_zc(h, signal_length);
	if(sil == 0  || signal_length <= 256) {
		free(h);
		return (num_ph - 1);
	}
	// if(signal_length < (glbl_window_width * 5)) {
	// 	return (num_ph - 1);
	// }
	int new_length = signal_length;
	// by width / div_interval ?
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
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;

	// double ste = short_time_energy(h, signal_length, glbl_window_width);
	float** complete_signal = (float**)malloc(sizeof(float*) * 3);
	float* signal = (float*)malloc(sizeof(float) * signal_length);
	for(int i = 0; i < signal_length; i++) {
		signal[i] = h[i];
	}
	
	int extra_size = floor(signal_length / glbl_window_width);

	complete_signal[1] = (float*)calloc(extra_size, sizeof(float));
	complete_signal[2] = (float*)calloc(extra_size, sizeof(float));
	for(int m = 0; m < extra_size; m++) {
		complete_signal[1][m] = f_cross_rate(&signal[m * glbl_window_width], glbl_window_width);
		complete_signal[2][m] = short_time_energy(&h[m * glbl_window_width], glbl_window_width, glbl_window_width); 
	}

	signal = mfcc(signal, signal_length, glbl_window_width, glbl_banks, glbl_paa_op);
	
	complete_signal[0] = signal;
	if(signal == NULL) {
		printf("Unable to produce mfcc for phoneme\n");
		for(int i = 1; i < j; i++) {
			phones[i]->score = DBL_MAX;
		}
		free(h);
		return - 1;
	}
	int new_size = mfcc_size(signal_length);
	if(new_size <= 0) {
		printf("No length found for a test sequence, possibly too short or window too large\nNo testing has been done\n length : %d :: || incr :: %d || mfcc_size :: %d\n", signal_length, glbl_window_width, new_size);
		for(int i = 1; i < j; i++) {
			phones[i]->score = DBL_MAX;
		}
		free(h);
		return -1;
	}
		
	
	DTW_ERROR = 0;
	best_so_far = DBL_MAX;
	int knn = 0;
	for(int i = 1; i < j; i++) {
		phones[i]->score = 0;
		// if(zc < ph_zc_max[i] && zc > ph_zc_min[i] && DTW_ERROR != NO_SIG_LEN &&
		//	ste < ste_max[i] && ste > ste_min[i]) {
		if(CLUST && !(KNN)) {
			dtw_clust(complete_signal, signal_length, phones[i], glbl_dtw_window);
		} else if(KNN && !(CLUST)) {
			knn = knn_mfccs_size(complete_signal, signal_length, glbl_k, NULL);
			break;
		} else if(KNN && CLUST) {
			knn = k_means(signal, signal_length, glbl_k);
			break;
		} else {
			dtw_frame(complete_signal, signal_length, phones[i], glbl_dtw_window);
		}
	}
	if(KNN) {
		for(int i = 1; i < num_ph; i++) {
			if(i == knn) {
				phones[i]->score = 0;
			} else {
				phones[i]->score = 1;
			}
		}
	}

	if(failed == 1) {
		printf("Failed Phoneme\n");
		failed = 0;
	} else {
		tested++;
	}
	double sml = DBL_MAX;
	int pos = 1;
	for(int i = 1; i < j; i++) {
		if(phones[i]->score < sml) {
			sml = phones[i]->score;
			pos = i;
		}
	}
	// export_results(phones[1]->index->name);
	free(h);
	free(complete_signal[0]);
	free(complete_signal[1]);
	free(complete_signal[2]);
	free(complete_signal);
	if(KNN) {
		if(phones[knn]->index->group_i == 0 && sil == 2) {
			return (num_ph - 1);
		}
		return knn;
	} else {
	     
		return pos;
	}
}

/** 
 * @brief Tests files knowing the boundary
 * 
 * @param h The signal to be classified
 * @param signal_length The length of @param h
 * @param p The previous phoneme, used for the grammar functions
 *
 * This is the main one-to-one testing function and can test phonemes using basic DTW, KNN, K-means, etc.
 * A 'complete_signal' can be produced which stores the STE and ZC values for each frame to be used in voice type KNN classification
 */
void test_phoneme(short* h, int signal_length, char* p)
{

	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	// int test_length = signal_length;

	int new_size = mfcc_size(signal_length);
	if(new_size <= 0) {
		printf("No length found for a test sequence, possibly too short or window too large\nNo testing has been done\n length : %d :: || incr :: %d || mfcc_size :: %d\n", signal_length, glbl_window_width, new_size);
		for(int i = 1; i < j; i++) {
			phones[i]->score = DBL_MAX;
		}
		return;
	}
	num_ph = j;

	// float ste = short_time_energy(h, signal_length, glbl_window_width);
	float** complete_signal = (float**)malloc(sizeof(float*) * 5);
	float* signal = (float*)malloc(sizeof(float) * signal_length);
	for(int i = 0; i < signal_length; i++) {
		signal[i] = h[i];
	}

	/* int sil = is_sil_ste(h, signal_length); */
	/* // int sil_zc = is_sil_zc(h, signal_length); */
	/* if(sil == 0  || signal_length <= 256) { */
	/* 	for(int i = 1; i < num_ph; i++) { */
	/* 		if(strcmp(p, phones[i]->index->name) == 0) { */
	/* 			sil_v_matrix[2][phones[i]->index->voice]++; */
	/* 		} */
	/* 	} */
	/* } */
	
	int start = 1, end = num_ph;
	
	/* int grp = knn_mfccs_voice_time(signal, test_length, glbl_voice_k, NULL); */
	/* if      (grp == 0) { start = 1;  end = 16; } */
	/* else if (grp == 1) { start = 16; end = 38; } */
	/* else               { start = 38; end = 39; } */

	/* int v_grp = 0; */
	/* for(int i = 1; i < num_ph; i++) { */
	/* 	if(strcmp(p, phones[i]->index->name) == 0) { */
	/* 		v_grp = phones[i]->index->voice; */
	/* 	} */
	/* } */
	
	/* sil_v_matrix[grp][v_grp]++; */
	/* printf("grp :: %d || v_grp :: %d\n", grp, v_grp); */
	/* int grp = knn_mfccs_group_time(signal, test_length, glbl_voice_k, NULL); */
	/* if      (grp == 0) { start = 1;  end = 7;       } */
	/* else if (grp == 1) { start = 7;  end = 9;       } */
	/* else if (grp == 2) { start = 9;  end = 16;      } */
	/* else if (grp == 3) { start = 16; end = 19;      } */
	/* else if (grp == 4) { start = 19; end = 24;      } */
	/* else if (grp == 5) { start = 24; end = 38;      } */
	/* else               { start = 38; end = 39;      } */
	/* int sil = is_sil_ste(h, signal_length); */
	/* int sil_zc = is_sil_zc(h, signal_length); */
	
	int extra_size = floor(signal_length / glbl_window_width);

	complete_signal[1] = (float*)calloc(extra_size, sizeof(float));
	complete_signal[2] = (float*)calloc(extra_size, sizeof(float));
	complete_signal[3] = (float*)calloc(extra_size, sizeof(float));
	complete_signal[4] = (float*)calloc(extra_size, sizeof(float));
	for(int m = 0; m < extra_size; m++) {
		complete_signal[1][m] = f_cross_rate(&signal[m * glbl_window_width], glbl_window_width);
		complete_signal[2][m] = short_time_energy(&h[m * glbl_window_width], glbl_window_width, glbl_window_width);
		complete_signal[3][m] = kurtosis(&signal[m * glbl_window_width], glbl_window_width);
		complete_signal[4][m] = flatness(&signal[m * glbl_window_width], glbl_window_width);
	}
	
	signal = mfcc(signal, signal_length, glbl_window_width, glbl_banks, glbl_paa_op);
	
	if(NORM) {
		normalise_mfcc(signal, new_size);
	}
	complete_signal[0] = signal;
	if(signal == NULL) {
		printf("Unable to produce mfcc for phoneme :: %s\n", p);
		for(int i = 1; i < j; i++) {
			phones[i]->score = DBL_MAX;
		}
		return;
	}
	
	DTW_ERROR = 0;
	best_so_far = DBL_MAX;
	int knn = 0;
	for(int i = start; i < end; i++) {
		phones[i]->score = 0;
		// if(stavg > stavg_noov[i] * 0.75 && stavg < stavg_noov[i] * 1.25) {
		if(DTW_ERROR != NO_SIG_LEN) { // &&
		   // ste < ste_max[i] && ste > ste_min[i]) {
			if(CLUST && !(KNN)) {
				dtw_clust(complete_signal, signal_length, phones[i], glbl_dtw_window);
			} else if(KNN && !(CLUST)) {
				knn = knn_mfccs_size(complete_signal, signal_length, glbl_k, p);
				break;
			} else if(KNN && CLUST) {
				knn = k_means(signal, signal_length, glbl_k);
				break;
			} else {
				dtw_frame(complete_signal, signal_length, phones[i], glbl_dtw_window);
			}
		} else {
			phones[i]->score = DBL_MAX;
		}
	}
	if(KNN || SVM) {
		for(int i = 1; i < num_ph; i++) {
			if(i == knn) {
				phones[i]->score = 0;
			} else {
				phones[i]->score = 1;
			}
		}
	}
	if(failed == 1) {
		printf("Failed Phoneme :: %s\n", p);
		failed = 0;
	} else {
		tested++;
	}
	export_results(p);
	free(complete_signal[0]);
	free(complete_signal[1]);
	free(complete_signal[2]);
	free(complete_signal);
		
	return;
}

/** 
 * @brief Exports the results of a test. Used for one-to-one testing, not boundary testing.
 * 
 * @param ph_code The code of the phoneme which has been tested.
 *
 * A result may be exported to correct, group or fail. If the classification is correct then 'correct' will be used, if it is was within the same
 * group then 'group' will be used, if neither then 'fail' will be used. This export will contain what the correct phoneme's score was and what the
 * classifications score was at the top with a list of all tested phonemes with their scores below.
 */
void export_results(char* ph_code)
{
	// init structs - load the proto sequence
	int j = 1, index = 0;
	double actual = 0, smallest = DBL_MAX;
	actual = actual;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;
	int t_index = 0;
	int zero_scores = 0;
	// short results[3] = {SHRT_MAX, SHRT_MAX, SHRT_MAX};
	// int s_indexes[3] = {0, 0, 0};
	/* for(int i = 1; i < j; i++) { */
	/* 	if(strcmp(phones[i]->index->name, ph_code) == 0) { */
	/* 		actual = phones[i]->score; */
	/* 		t_index = i; */
	/* 		if(phones[i]->score == 0) { */
	/* 			zero_scores++; */
	/* 		} */
	/* 	} */
	/* 	if(fabs(phones[i]->score) < fabs(smallest)) { */
	/* 		smallest = phones[i]->score; */
	/* 		index = i; */
	/* 	}	 */
	/* } */
	for(int i = 1; i < j; i++) {
		if(strcmp(phones[i]->index->name, ph_code) == 0) {
			actual = phones[i]->score;
			t_index = i;
			/* if(phones[i]->score == 0) { */
			/* 	smallest = phones[i]->score; */
			/* 	index = i; */
			/* 	sil_corr++; */
			/* 	break; */
			/* } */
		}
		if(phones[i]->score == 0) { 
			zero_scores++; 
		}
		if(fabs(phones[i]->score) < fabs(smallest)) {
			smallest = phones[i]->score;
			index = i;
		}	
	}

	char file_name[25] = {0};
	if(index == 0) {
		index = num_ph - 1;
	}
	if(strcmp(ph_code, "h#") == 0
	   && (strcmp(phones[index]->index->name, "epi") == 0 || strcmp(phones[index]->index->name, "pau") == 0)) {
		correct++;
		strcat(file_name, "./correct/");
	} else if(strcmp(ph_code, "epi") == 0
		  && (strcmp(phones[index]->index->name, "h#") == 0 || strcmp(phones[index]->index->name, "pau") == 0)) {
		correct++;
		strcat(file_name, "./correct/");
	} else if(strcmp(ph_code, "pau") == 0
		  && (strcmp(phones[index]->index->name, "epi") == 0 || strcmp(phones[index]->index->name, "h#") == 0)) {
		correct++;
		strcat(file_name, "./correct/");
	} else if(zero_scores > 1) {
		zero_fails++;
		fails++;
		strcat(file_name, "./fails/");	
	} else if(strcmp(ph_code, phones[index]->index->name) == 0) {
		correct++;
		strcat(file_name, "./correct/");
	} else if(strcmp(phones[index]->index->group, phones[t_index]->index->group) == 0) {
		group++;
		strcat(file_name, "./group/");
	} else {
		fails++;
		strcat(file_name, "./fails/");
	}
	strcat(file_name, ph_code);
	strcat(file_name, ".txt");
	FILE* fp = fopen(file_name, "a");
	// printf("test filename = %s\n", file_name);
	if(fp == NULL) {
		printf("Failed saving phoneme result\n");
		exit(-1);
	}
	// printf("index = %d\n", s_indexes[0]);
	// printf("results = %d\n", results[0]);
	fprintf(fp, "==== testing %s (pos %d, group %s) and got %s (pos %d, group %s) with a score of %f | actual phone got a score of %f \n", ph_code, t_index, phones[t_index]->index->group, phones[index]->index->name, phones[index]->index->i, phones[index]->index->group, smallest, actual);
	for(int i =1; i < j; i++) {
		if(phones[i]->score != DBL_MAX) {
			fprintf(fp, "Score : %f | ph_code : %s\n", phones[i]->score, phones[i]->index->name);
		}
		phones[i]->score = 0;
	}
	// printf("t_index %d | index %d | num_ph %d\n", t_index, index, num_ph);
	matrix[t_index][index]++;
	memset(file_name, 0, sizeof(file_name));
	fclose(fp);
	for(int i = 1; i < j; i++) {
		phones[i]->score = 0;
	}
	done++;
	if(done % (BOUNDS ? 500 : 2500) == 0)
		printf("::     COMPLETED TEST     ::  %02d:%02d:%02d  ::  %05d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), done);
	return;
}

/** 
 * @brief Tests a found phoneme using boundary detection.
 * 
 * @param filename The filename of the reference file
 * @param sequence The found phoneme audio sequence
 * @param length The length of @param sequence
 * @param offset The offset of the file, used when NOSIL is used as the start of testing is no longer the start of the .wav file due to the removal of the
 * leading and trailing silence.
 */
void utterance_test(char* filename, short* sequence, int length, int offset)
{
	char file[128];
	if(MALE) {
		memcpy(file, male_dir, sizeof(char) * (strlen(male_dir) + 1));
	} else if (FEMALE) {
		memcpy(file, female_dir, sizeof(char) * (strlen(female_dir) + 1));
	} else if (SPKR1 && !(CAR || STRT || CAFE || WHITE)) {
		memcpy(file, spkr1_dir, sizeof(char) * (strlen(spkr1_dir) + 1));
	} else if (SPKR1 && CAR) {
		memcpy(file, spkr1_car_dir, sizeof(char) * (strlen(spkr1_car_dir) + 1));
	} else if (SPKR1 && STRT) {
		memcpy(file, spkr1_strt_dir, sizeof(char) * (strlen(spkr1_strt_dir) + 1));
	} else if (SPKR1 && CAFE) {
		memcpy(file, spkr1_cafe_dir, sizeof(char) * (strlen(spkr1_cafe_dir) + 1));
	} else if (SPKR1 && WHITE) {
		memcpy(file, spkr1_white_dir, sizeof(char) * (strlen(spkr1_white_dir) + 1));
	} else if (SPKR1_NOSIL) {
		memcpy(file, spkr1_nosil_dir, sizeof(char) * (strlen(spkr1_nosil_dir) + 1));
	} else {
		memcpy(file, test_dir, sizeof(char) * (strlen(test_dir) + 1));
	}
	strcat(file, "res/");
	strncat(file, filename, strlen(filename) - 4);
	strcat(file, ".res");
	FILE* fp = fopen(file, "a");
	if(!fp) { printf("Error opening res '%s' file\n", file); }
	int bound = 0, res = 0;
	int current_pos = 0;
	while(1) {
		bound = next_boundary(sequence, length);
		if(bound == 0) {
			sequence = shift_and_reduce(sequence, length, 32);
			length = (length - 32);
			continue;
		}
		if(bound == -1) {
			break;
		}
		fprintf(fp, "%d ", current_pos + offset);
		current_pos += bound;
		fprintf(fp, "%d ", current_pos + offset);
		short* h = calloc(bound, sizeof(short));
		for(int i = 0; i < bound; i++) {
			h[i] = sequence[i];
		}
		res = test_phoneme_utterance(h, bound);
		if(res == -1) {
			printf("none found...\n");
		} else {
			fprintf(fp, "%s \n", phones[res]->index->name);
		}
		sequence = shift_and_reduce(sequence, length, bound);
		length = (length - bound);
	}
	free(sequence);
	fclose(fp);
}

/** 
 * @brief Merges any contiguous silence classifications into one
 * 
 * @param filename The name of the file to find and merge silences
 *
 * Any classifications of subsequent silences results in the same outcome for the user but a different WER and so it is decided that as one or more silences
 * acts the same then they can be combined together.
 */
void merge_silences(char* filename)
{
	char line[256];
	char* end_ptr;
	char *p; 
	int amount = 0;
	FILE* fp = fopen(filename, "r");
	if(!fp) {
		printf("Error opening file :: %s when merging silences\n", filename);
		return;
	}
	while(fgets(line, sizeof(line), fp)) {
		amount++;
	}
	struct guess* gs = (struct guess*)malloc(amount * sizeof(struct guess));
	int n = 0;
	fseek(fp, 0, SEEK_SET);
	while(fgets(line, sizeof(line), fp)) {
		p = strtok(line," ");
		gs[n].start = strtol(p, &end_ptr, 10);
	      	p = strtok(NULL, " ");
		gs[n].end = strtol(p, &end_ptr, 10);
		p = strtok(NULL, " ");
		gs[n].name = (char*)malloc((strlen(p) + 1) * sizeof(char));
		memcpy(gs[n].name, p, (strlen(p) + 1) * sizeof(char));
		n++;
	}

	for(int i = 1; i < amount; i++) {
		if( strcmp(gs[i-1].name, "sil") == 0 && strcmp(gs[i].name, "sil") == 0) {
			free(gs[i-1].name);
			gs[i-1].name = NULL;
			gs[i].start = gs[i-1].start;
		}
	}
	fp = freopen(filename, "w", fp);
	for(int i = 0; i < amount; i++) {
		if(gs[i].name == NULL)
			continue;
		fprintf(fp, "%ld %ld %s\n", gs[i].start, gs[i].end, gs[i].name);
	}
	fclose(fp);
	for(int i = 0; i < amount; i++) {
		if(gs[i].name != NULL) {
			free(gs[i].name);
		}
	}
	free(gs);

	return;
}

/** 
 * @brief Calculates the size of all MFCCs currently in use
 * 
 * 
 * @return The current size of all MFCCs
 */
float data_size(void)
{
	float bytes_saved = 0;
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;
	}
	num_ph = j;
	for(int i = 1; i < num_ph; i++) {
		for(int m = 0; m < phones[i]->size_count; m++) {
			if(phones[i]->size[m] == 0)
				continue;
			bytes_saved += (phones[i]->size[m] * sizeof(float));
		}
	}
	
	return (bytes_saved / 1024 / 1024);
}

/** 
 * Removes any sil MFCCs above a certain frame length as they do not vary much among frames, removes coefficients with minimal variance among phonemes
 * and reduces coefficients with low, but not minimal, variance to a short rather than a float to save space.
 * 
 * 
 * @return The amount of data reduced by the function
 */
float reduce_data(void)
{
	float bytes_saved = 0;
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;
	}
	num_ph = j;
	int trunc = floor(glbl_banks * glbl_test_trunc);
	shorted = calloc((trunc + 1), sizeof(int));
	removed = calloc((trunc + 1), sizeof(int));
	
	for(int t = 0; t < trunc; t++) {
		long double variance = 0, mean = 0, sum = 0, n = 0, std = 0;
		for(int i = 1; i < num_ph - 3; i++) {
			for(int m = 0; m < phones[i]->size_count; m++) {
				if(phones[i]->size[m] == 0)
					continue;
				int co = 0;
				for(int p = 0; p < phones[i]->size[m]; p++) {
					if(co == t) {
						sum += phones[i]->mfcc[m][p];
						n++;
					}
					co++;
					if(co >= trunc) {
						co = 0;
					}
				}	
			}
		}
		mean = (sum / n);
		for(int i = 1; i < num_ph - 3; i++) {
			for(int m = 0; m < phones[i]->size_count; m++) {
				if(phones[i]->size[m] == 0)
					continue;
				int co = 0;
				for(int p = 0; p < phones[i]->size[m]; p++) {
					if(co == t) {
						variance += pow( (phones[i]->mfcc[m][p] - mean), 2);
						std += pow( (phones[i]->mfcc[m][p] - mean), 2);
					}
					co++;
					if(co >= trunc) {
						co = 0;
					}
				}	
			}
		}
		variance /= n;
		std /= n;
		std = sqrt(std);
		printf("Variance for coefficient %d is %Lf ", t, variance);
		printf("Standard Deviation for coefficient %d is %Lf ", t, std);
		printf("Mean for coefficient is %Lf\n", mean);
	}
	
	for(int i = 1; i < num_ph; i++) {
		for(int m = 0; m < phones[i]->size_count; m++) {
			if(phones[i]->size[m] == 0)
				continue;
			if(phones[i]->used[m] <= 1) {
				bytes_saved += (phones[i]->size[m] * sizeof(float));
				free(phones[i]->mfcc[m]);
				phones[i]->mfcc[m] = NULL;
				phones[i]->size[m] = 0;
				phones[i]->reduced_count--;
			} else if(strcmp(phones[i]->index->name, "sil")  == 0) {
				if(phones[i]->size[m] > 1024) {
					bytes_saved += (phones[i]->size[m] * sizeof(float));
					free(phones[i]->mfcc[m]);
					phones[i]->mfcc[m] = NULL;
					phones[i]->size[m] = 0;
					phones[i]->reduced_count--;
				}
			} else {
				int co = 0;
				for(int p = 0; p < phones[i]->size[m]; p++) {
					if(co > 18) {
						phones[i]->mfcc[m][p] = (signed char)phones[i]->mfcc[m][p];
						bytes_saved += (sizeof(float) - sizeof(signed char));
					}
					// co++;
					if(co >= trunc) {
						co = 0;
					}
				}
			}
		}
	}

	return (bytes_saved / 1024 / 1024);
}

/** 
 * @brief Outputs the boundary test results to the screen, and the results with which paramters to a file
 * 
 * @param m_fp The file to output the results to.
 */
void bounds_output(FILE* m_fp)
{
	// float avg_len = (total_lengths / tested_files);
	// float avg_dist = (min_edit / worst) * 100;
	long double end_wer = ((ins+subs+delts)/(refs));

	// printf("::        MIN_DIST        ::  %02d:%02d:%02d  ::  %.2f%%:%.2f\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), avg_dist, avg_len);
	printf("::     COR/INS/SUB/DEL    ::  %02d:%02d:%02d  ::  %.0Lf:%.0Lf:%.0Lf:%.0Lf\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), corr, ins, subs, delts);
	printf("::        AVERAGES        ::  %02d:%02d:%02d  ::  %.2Lf:%.2Lf:%.2Lf:%.2Lf\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), (corr/total_dists) * 100, (ins/total_dists) * 100, (subs/total_dists) * 100, (delts/total_dists) * 100);
	printf("::           WER          ::  %02d:%02d:%02d  ::  %.2Lf\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), end_wer);
	printf("::           ...          ::  %02d:%02d:%02d  ::  %.2Lf:%.2Lf:%.2Lf\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), avg_wer/tested_files, low_wer, high_wer);
	printf("::        BEST FILE       ::  %02d:%02d:%02d  ::  %s : %.2Lf\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), &best_file[12], low_wer);
	printf("::       WORST FILE       ::  %02d:%02d:%02d  ::  %s : %.2Lf\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), &worst_file[12], high_wer);
	// printf("::     TOT/DEL/INS/SUB    ::  %02d:%02d:%02d  ::  %Lf:%Lf:%Lf:%Lf\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), total_lengths, delts, ins, subs);

	// float avg_len = (total_lengths / tested_files) * 100;
	float avg_dist = (min_edit / total_dists) * 100;
	fprintf(m_fp, "time = %ld \n", time(NULL));
	fprintf(m_fp, "PAA_OP %d || PAA :: %d || BANKS :: %d || WINDOWS :: %d || LIMIT :: %d || INTERVAL :: %d || NFFT : %d || MFCCS : %d || D :: %d || DD :: %d || TRUNC :: %.2f || EXTRA :: %d || CLUST :: %d || ITER :: %d\n\n", glbl_paa_op, glbl_paa, glbl_banks, glbl_window_width, glbl_dtw_window, glbl_interval_div, glbl_nfft, glbl_mfcc_num, DELTA, DELTA_DELTA, floor(glbl_banks* glbl_test_trunc), EXTRA, glbl_clust_num, glbl_test_iter);
	fprintf(m_fp, "window_width = %d \n banks = %d \n paa = %d\n dtw_window = %d \n nfft = %d \n interval div = %d  \n limits changed = %d\n", glbl_window_width, glbl_banks, glbl_paa, glbl_dtw_window, glbl_nfft, glbl_interval_div, limit_changed);
	fprintf(m_fp, "zc_incr = %d \n ste_incr = %d \n entr_incr = %d \n larg_incr = %d \n neg_incr = %d \n", glbl_zc_incr, glbl_ste_incr, glbl_entr_incr, glbl_larg_incr, glbl_neg_incr); 

	fprintf(m_fp, "NORMAL:: \ncorrect = %Lf \n insertions = %Lf \n deletions = %Lf \n substitutions = %Lf \n percentage : %.2f%% \n WER : %.2Lf\n", corr, ins, delts, subs, avg_dist, ((ins+subs+delts)/refs) * 100);

	FILE* percents = fopen("../Dynamic_Time_Warping/percentages.txt", "a");
	fprintf(percents, "%.2Lf : %.2Lf : %.2Lf : %.2Lf : %.2Lf : %d : %d : %d : %d : %d\n", end_wer, corr, ins, subs, delts,
		glbl_zc_incr, glbl_ste_incr, glbl_entr_incr, glbl_larg_entr_incr, glbl_larg_ste_incr);
	
	fclose(percents);

	return;
}

/** 
 * @brief Outputs the confusion matrix and percentage correct files for use with the given Python scripts.
 * 
 * @param m_fp The file to output the matrix and percentages to.
 *
 * The percentages correct are defined as the amount correct over the amount observed.
 */
void matrix_output(FILE* m_fp)
{

	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;
	}
	num_ph = j;

	for(int i = 1; i < num_ph; i++) {
		fprintf(m_fp, "%s | ", p_codes[i]);
		for(int j = 1; j < num_ph; j++) {
			fprintf(m_fp," %d |", matrix[i][j]);
		}
		fprintf(m_fp, "\n");
	}

	// phoneme matrix for python
	FILE* py = fopen("../Dynamic_Time_Warping/confusion_matrix.txt", "w");
	for(int i = 1; i < num_ph; i++) {
		fprintf(py, "%s | ", p_codes[i]);
		for(int j = 1; j < num_ph; j++) {
			fprintf(py," %d |", matrix[i][j]);
		}
		fprintf(py, "\n");
	}
	fclose(py);

	FILE* perc = fopen("../Dynamic_Time_Warping/per_correct.txt", "w");
	for(int i = 1; i < num_ph; i++) {
		if(per_correct[i] != 0) {
			fprintf(perc, "%s | %.2f\n", p_codes[i], ((float)matrix[i][i] / (float)per_correct[i] * 100));
		}
	}
	fclose(perc);
	
	FILE* percents = fopen("../Dynamic_Time_Warping/percentages.txt", "a");
	float per = ((float)correct / ((float)group + (float)correct + (float)fails)) * 100.0;
	fprintf(percents, "%.2f : %d\n", per, glbl_dtw_window);
	
	fclose(percents);
	// group matrix
	for(int i = 0; i < 7; i++) {
		fprintf(m_fp, "%s | ", p_group[i]);
		for(int j = 0; j < 7; j++) {
			fprintf(m_fp," %d |", group_matrix[i][j]);
		}
		fprintf(m_fp, "\n");
	}
	// voice matrix
	char* voices[] = { "vl", "vd", "nv" };
	for(int i = 0; i < 3; i++) {
		fprintf(m_fp, "%s | ", voices[i]);
		for(int j = 0; j < 3; j++) {
			fprintf(m_fp," %d |", voice_matrix[i][j]);
		}
		fprintf(m_fp, "\n");
	}

	for(int i = 0; i < 3; i++) {
		fprintf(m_fp, "%s | ", voices[i]);
		for(int j = 0; j < 3; j++) {
			fprintf(m_fp," %d |", sil_v_matrix[i][j]);
		}
		fprintf(m_fp, "\n");
	}

	fprintf(m_fp, "\n");

	return;
}

/** 
 * @brief Outputs the one-to-one testing to the user.
 * 
 * @param m_fp The file to output the results to.
 */
void test_output(FILE* m_fp)
{
	float per = 0;
	per = ((float)correct / ((float)group + (float)correct + (float)fails)) * 100.0;
	fprintf(m_fp, "time = %ld \n", time(NULL));
	fprintf(m_fp, "PAA_OP %d || PAA :: %d || BANKS :: %d || WINDOWS :: %d || LIMIT :: %d || INTERVAL :: %d || NFFT : %d || MFCCS : %d || D :: %d || DD :: %d || TRUNC :: %.2f || EXTRA :: %d || K :: %d\n\n", glbl_paa_op, glbl_paa, glbl_banks, glbl_window_width, glbl_dtw_window, glbl_interval_div, glbl_nfft, glbl_mfcc_num, DELTA, DELTA_DELTA, glbl_test_trunc, EXTRA, glbl_k);
	fprintf(m_fp, "window_width = %d \n banks = %d \n paa = %d\n dtw_window = %d \n nfft = %d \n interval div = %d  \n limits changed = %d\n", glbl_window_width, glbl_banks, glbl_paa, glbl_dtw_window, glbl_nfft, glbl_interval_div, limit_changed);

	fprintf(m_fp, "NORMAL:: \ncorrect = %d \n fails = %d \n groups = %d\n sil_correct = %d\n percentage : %.2f%%\n\n", correct, fails, group, sil_corr, per);

	FILE* percents = fopen("../Dynamic_Time_Warping/percentages.txt", "a");
	per = ((float)correct / ((float)group + (float)correct + (float)fails)) * 100.0;
	fprintf(percents, "%.2f : %d : %f\n", per, glbl_banks, floor(glbl_test_trunc * glbl_banks));
	
	fclose(percents);
	
	return;
}

/** 
 * @brief Trims the silence data from a given .wav file
 * 
 * @param fp The reference file associated with the .wav file
 * @param sequence The associated .wav file in an array
 * 
 * @return The new wav structure with silence removed, the size updated and the offset
 *
 * When using a NOSIL dataset this function returns the structure with the silence removed and provides the offset so that the boundary .res files are correct.
 */
struct wav_file* trim_silence(FILE* fp, short* sequence)
{
	fseek(fp, 0, SEEK_SET);
	long start = 0, end = 0, f_start = 0, f_end = 0, length = 0;
	char line[256];
	char* end_ptr;
	char *p;
	int n = 0;
	while(fgets(line, sizeof(line), fp)) {
		p = strtok(line," ");
		start = strtol(p, &end_ptr, 10);
	      	p = strtok(NULL, " ");
		end = strtol(p, &end_ptr, 10);
		p = strtok(NULL, " ");
		while(p[strlen(p) - 1] == ' ' || p[strlen(p) - 1] == '\n') {
			p[strlen(p) - 1] = '\0';
		}
		if(n == 0) {
			f_start = start;
		}
		n++;
	}
	f_end = end;
	length = (f_end - f_start);
	struct wav_file* result = (struct wav_file*)malloc(sizeof(struct wav_file));
	result->seq = (short*)calloc(length, sizeof(short));
	result->length = length; 
	for(int i = f_start; i < f_end; i++) {
		result->seq[i - f_start] = sequence[i];
	}
	result->offset = f_start;
	free(sequence);
	fseek(fp, 0, SEEK_SET);
	return result;
	
}

/** 
 * @brief The main control function of the testing process.
 *
 * Determines and sets the test dataset file path, iterates of the files in the dataset file path, produces and outputs the results of the testing process.
 */
void test(void)
{
	reset();
	done = 0;
	if(BOUNDS) {
		best_file = (char*)calloc(1, sizeof(char));
		worst_file = (char*)calloc(1, sizeof(char));
	}
	DIR *p;
	struct dirent *pp;
	char dir_name[35];
	char res_dir_name[55];
	if(MALE) {
		memcpy(dir_name, male_dir, sizeof(char) * (strlen(male_dir) + 1));
		memcpy(res_dir_name, res_male_dir, sizeof(char) * (strlen(res_male_dir) + 1));
	} else if (FEMALE) {
		memcpy(dir_name, female_dir, sizeof(char) * (strlen(female_dir) + 1));
		memcpy(res_dir_name, res_female_dir, sizeof(char) * (strlen(res_female_dir) + 1));
	} else if (SPKR1 && !(CAR || STRT || CAFE || WHITE)) {
		memcpy(dir_name, spkr1_dir, sizeof(char) * (strlen(spkr1_dir) + 1));
		memcpy(res_dir_name, res_spkr1_dir, sizeof(char) * (strlen(res_spkr1_dir) + 1));
	} else if (SPKR1 && CAR) {
		memcpy(dir_name, spkr1_car_dir, sizeof(char) * (strlen(spkr1_car_dir) + 1));
		memcpy(res_dir_name, res_spkr1_car_dir, sizeof(char) * (strlen(res_spkr1_car_dir) + 1));
	} else if (SPKR1 && STRT) {
		memcpy(dir_name, spkr1_strt_dir, sizeof(char) * (strlen(spkr1_strt_dir) + 1));
		memcpy(res_dir_name, res_spkr1_strt_dir, sizeof(char) * (strlen(res_spkr1_strt_dir) + 1));
	} else if (SPKR1 && CAFE) {
		memcpy(dir_name, spkr1_cafe_dir, sizeof(char) * (strlen(spkr1_cafe_dir) + 1));
		memcpy(res_dir_name, res_spkr1_cafe_dir, sizeof(char) * (strlen(res_spkr1_cafe_dir) + 1));
	} else if (SPKR1 && WHITE) {
		memcpy(dir_name, spkr1_white_dir, sizeof(char) * (strlen(spkr1_white_dir) + 1));
		memcpy(res_dir_name, res_spkr1_white_dir, sizeof(char) * (strlen(res_spkr1_white_dir) + 1));
	} else if (SPKR1_NOSIL) {
		memcpy(dir_name, spkr1_nosil_dir, sizeof(char) * (strlen(spkr1_nosil_dir) + 1));
		memcpy(res_dir_name, res_spkr1_nosil_dir, sizeof(char) * (strlen(res_spkr1_nosil_dir) + 1));
	} else {
		memcpy(dir_name, test_dir, sizeof(char) * (strlen(test_dir) + 1));
		memcpy(res_dir_name, res_test_dir, sizeof(char) * (strlen(res_test_dir) + 1));
	}
	if(!(p = opendir (dir_name))) {
		printf("Failed to open folder %s\n", dir_name);
		return;
	}
	char pathname[1024];
	FILE* fp;
	
	int files_to_test = 0;
	if(p != NULL) {
		while((pp = readdir(p)) != NULL) {
			files_to_test++;
		}
	}
	(void) closedir (p);
	if(!(p = opendir (dir_name))) {
		printf("Failed to open folder %s after counting files\n", dir_name);
		return;
	}
	
	int chunk = files_to_test;
	if(SPLIT_DATA) {
		chunk = floor(files_to_test / glbl_test_iter);
		files_to_test = chunk;
	}
	int current_file = 0;
	printf(":: Testing folder :: %s :: %d files to test\n", dir_name, files_to_test);
	if(p != NULL) {
		while ((pp = readdir (p)) != NULL) {
			if(current_file < current_chunk) {
				current_file++;
				continue;
			}
			if(current_file >= (current_chunk + chunk)) {
				current_chunk += chunk;
				break;
			}
			current_file++;
			// printf("current file: %d :: current_chunk : %d :: chunk :: %d\n", current_file, current_chunk, chunk);
			short* sequence  = NULL;
			int length = strlen(pp->d_name);
			struct wav_file* wv = NULL;
			int seq_len = 0;
			if (strncmp(pp->d_name + length - 4, ".wav", 4) == 0 || strncmp(pp->d_name + length - 4, ".WAV", 4) == 0) {
				int len = 4;
				sprintf(pathname, "%s%s", dir_name, pp->d_name);
				fp = fopen(pathname, "rb");
				if(!fp) { printf("Error opening wav '%s' file\n", pathname); continue; }
				wv = read_wav(fp);
				sequence = wv->seq;
				seq_len = wv->length;
			        pp->d_name[length-len] = '\0';
				strcat(pp->d_name, ".PHN");
				sprintf(pathname, "%s%s", dir_name, pp->d_name);
				fp = fopen(pathname, "r");
				if(!fp) { printf("Error opening phone '%s' file\n", pathname); continue; }
				if(BOUNDS) {
					int offset = 0;
					struct wav_file* trimmed = NULL;
					if(SPKR1_NOSIL) {
						trimmed = trim_silence(fp, sequence);
						sequence = trimmed->seq;
						offset = trimmed->offset;
						seq_len = trimmed->length;
					}
					utterance_test(pp->d_name, sequence, seq_len, offset);
					char res_pathname[1024];
					pp->d_name[length-len] = '\0';
					strcat(pp->d_name, ".res");
					sprintf(res_pathname, "%s%s", res_dir_name, pp->d_name);
					merge_silences(res_pathname);
					minimum_edit_distance(res_pathname, pathname);
					tested_files++;
					if(trimmed != NULL) {
						free(trimmed);
					}
					if(tested_files % 250 == 0)
						printf("::     COMPLETED TEST     ::  %02d:%02d:%02d  ::  %05d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), tested_files);
				} else {
					allocate_ph(fp, sequence, 1);
					free(sequence);
				}
			}
			free(wv);
			
		}
		(void) closedir (p);
	}

	FILE* m_fp = fopen("../Testing/TEST/matrix.txt", "a");
	if(m_fp == NULL) {
		printf("Failed saving phoneme result\n");
		exit(-1);
	}
	
	if(BOUNDS) {
		bounds_output(m_fp);
		char best_file_out[512];
		char best_file_in[512];
		int slash_count = 0;
		sprintf(best_file_out, "%s-%d_%d_%d_%d_%d-%.2Lf.txt", best_file, glbl_zc_incr, glbl_ste_incr,
			glbl_entr_incr, glbl_neg_incr, glbl_larg_incr, low_wer);
		while(slash_count != 4) {
			best_file++;
			if(best_file[0] == '/')
				slash_count++;
		}
		sprintf(best_file_out, "../Testing/Best_Bounds/%s-%d_%d_%d_%d_%d.txt", best_file, glbl_zc_incr, glbl_ste_incr,
			glbl_entr_incr, glbl_neg_incr, glbl_larg_incr);
		
		FILE* best_out = fopen(best_file_out, "w");
		if(best_out == NULL) {
			printf("Failed to open out file %s\n", best_file_out);
			exit(-1);
		}
		FILE* best_in = fopen(best_file_in, "r");
		if(best_in == NULL) {
			printf("Failed to open out file %s\n", best_file_in);
			exit(-1);
		}
		char a;
		while ( (a = fgetc(best_in)) != EOF ) {
			fputc(a, best_out);
		}
		fclose(best_out);
		fclose(best_in);
	}


	if(!BOUNDS) {
		test_output(m_fp);
		matrix_output(m_fp);
	}

	fclose(m_fp);

	FILE* zc_fp = fopen("../Testing/TEST/zero_cross.txt", "a");
	if(zc_fp == NULL) {
		printf("Failed saving zero cross result\n");
		exit(-1);
	}
	// float avg = 0;
	for(int p = 1; p < num_ph; p++) {
		// avg = (ph_zc[p] / zc_added[p]);
		// fprintf(zc_fp, "Phoneme :: %s || Total ZC :: %d || Count :: %d || Avg :: %.3f || Min :: %d || Max :: %d || ST Avg :: %f || ST Avg NoOv :: %f\n", p_codes[p], ph_zc[p], zc_added[p], avg, ph_zc_min[p], ph_zc_max[p], zc_st_avg[p], zc_st_avg_no_oc[p]);
		fprintf(zc_fp, "Phoneme :: %s || ST Min :: %f || ST Max :: %f\n", p_codes[p], ph_zc_min[p], ph_zc_max[p]);
	}
	fprintf(zc_fp, "\n");

	for(int p = 1; p < num_ph; p++) {
		// avg = (ph_zc[p] / zc_added[p]);
		// fprintf(zc_fp, "Phoneme :: %s || Total ZC :: %d || Count :: %d || Avg :: %.3f || Min :: %d || Max :: %d || ST Avg :: %f || ST Avg NoOv :: %f\n", p_codes[p], ph_zc[p], zc_added[p], avg, ph_zc_min[p], ph_zc_max[p], zc_st_avg[p], zc_st_avg_no_oc[p]);
		fprintf(zc_fp, "Phoneme :: %s || STE Min :: %f || STE Max :: %f\n", p_codes[p], ste_min[p], ste_max[p]);
	}
	
	fclose(zc_fp);




	
	

	return;
}

