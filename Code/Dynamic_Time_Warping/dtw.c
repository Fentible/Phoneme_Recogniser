/**
 * @file   dtw.c
 * @author T. Buckingham
 * @date   Wed May  4 15:41:32 2022
 * 
 * @brief  The main control file of the program which also contains the DTW functions
 * 
 * This file contains program dependant configuration variables along with
 * all required dynamic time warping files. This file was decided to contain these and \fn main()
 * due to the project originally focusing exclusively on DTW as a classfification method. 
 */

#define _GNU_SOURCE
#define _USE_MATH_DEFINES

#include "dtw.h"

int glbl_banks = 40;             /* The number of filter banks used for the MFCC */
int glbl_window_width = 128;     /* The width of the Hanning window */
int glbl_paa_op = 1;             /* Whether PAA is used; 1 for no, 0 for yes */
int glbl_paa = 2;                /* The divison of decimation for PAA */
int glbl_dtw_window = 200;       /* The DTW window limit as x / 1000 to produce a percentage */
int glbl_interval_div = 2;       /* The overlap division for \var glbl_window_width */
int glbl_nfft = 512;             /* The NFFT used when creating the filter bank */
int glbl_mfcc_num = 25;          /* The MFCC limit for each phoneme */
int glbl_clust_num = 0;          /* The number of cluster centroids */
int glbl_k = 7;                  /* The k value for base KNN */
float glbl_test_trunc = 24;      /* The number of coefficients kept from the filter bank */
int no_window = 0;              
int glbl_test_iter = 1;          /* The number of testing iterations to perform */
int glbl_group_k = 7;            /* The k value for group KNN */
int glbl_voice_k = 7;            /* The k value for voiced KNN */
int glbl_frame_limit = INT_MAX;  /* The MFCC frame limit */

int glbl_zc_incr;                /* The zero cross threshold amount as an absolute difference */
int glbl_ste_incr;               /* The short time energy threshold as a percentage difference */
int glbl_entr_incr;              /* The entropy threshold as a percentage difference */
int glbl_neg_incr;               /* The negative threshold for entropy (DEPRECATED) */
int glbl_larg_incr;              /* The large entropy threshold (DEPRECATED) */
int glbl_larg_entr_incr;         /* A larger entropy threshold used for boundary detection as a percentage */
int glbl_larg_ste_incr;          /* A larger short time energy threshold used for boundary detection as a percentage */

int largest_mfcc = 0;            /* Stores the largest MFCC in frames. Used for output */
int largest_index = 0;           
float largest_value = 0;         /* The largest value within the MFCCs. Used for output */
float smallest_value = FLT_MAX;  /* The smallest value within the MFCCs. Used for output */

short failed = 0;                /* Used to check if a DTW functions returns with an error */

int DTW_ERROR;                   /* If there was an error performing DTW */
int DELTA       = 0;             /* If delta coefficients should be used */
int DELTA_DELTA = 0;             /* If delta-delta coefficients should be used */
int NORM        = 0;             /* If the features should be normalised using min-max */
int SIMPLE      = 0;             /* If the simple, normal, version of MFCCs should be used */
int EXTRA       = 0;             /* If ZC and STE should be used during classification */
int CLUST       = 0;             /* If clustering should be performed during traning and testing */
int BOUNDS      = 0;             /* If testing should use boundary detection */
int LOG_E       = 0;             /* If additional information such as the log entropy should be added to the end of the MFCCs */
int AVG         = 0;             /* If the MFCCs should be averaged */
int KNN         = 0;             /* If KNN should be used in classification */
int GROUP       = 0;             /* If group KNN should be used in classification */
int VOICED      = 0;             /* If voice KNN should be used in classsifcation */
int GRAM        = 0;             /* If the grammar methods should be used in training and testing */
int PCA         = 0;             /* If the MFCCs should be exported to be used in PCA Python scripts; no testing will be done */
int ZC          = 0;             /* If frame-by-frame zero cross should be generated and used for voiced KNN */
int STE         = 0;             /* If frame-by-frame short time energy should be generated and used for voiced KNN */
int MEAN_SIZE   = 0;             /* If the MFCCs should be averaged by their length in frames */
int KURT        = 0;             /* If kurtosis should be used for group KNN */
int ENTR        = 0;             /* If entropy should be used for group KNN */
int SPKR1       = 0;             /* If Speaker 1's data should be used for training and testing */
int SPKR1_NOSIL = 0;             /* If Speaker 1 (No silence) data should be used for training and testing */
int WHITE       = 0;             /* If the white noise version of Speaker 1's data should be used for testing */
int STRT        = 0;             /* If the street noise version of Speakers 1's data should be used for testing */
int CAR         = 0;             /* If the car noise version of Speakers 1's data should be used for testing */
int CAFE        = 0;             /* If the cafe noise version of Speakers 1's data should be used for testing */
int SPLIT_DATA  = 0;             /* If the testing data should be divided by the number of testing iterations */
int SVM         = 0;             /* If SVM should be used for classification; the SVM weights should already be produced */
int EXPORT      = 0;             /* Exports the phonemes for use on the device */
int THREAD      = 0;             /* Uses threads to create the MFCCs */
int RAW         = 0;             /* Uses raw time domain signals for some classifications */

int MALE   = 0;                  /* If the male data should be used for testing */
int FEMALE = 0;                  /* If the female data should be used for testing */

int trained   = 0;               /* The number of raw audio phoneme signals added to the model */
int tested    = 0;               /* The number of tests performed */
int mfcced    = 0;               /* The number of MFCCs creared */
int clustered = 0;               /* The number of clusters created */
int exported  = 0;               /* The number of MFCCs exported */

double best_so_far = DBL_MAX;    /* The best so far score for DTW; not used during KNN */

time_t avg_test_time          = 0; /* The average test time per phoneme */
long double total_test_time   = 0; /* The sum of all test times during testing */
long double total_dtw_tests   = 0; /* The total number of tests performed during testing */

time_t start = 0;                  /* The program start time; used for output */

pthread_t* threads = NULL;         
int thread_count = 0;
static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

int clusters_left = 0;

int max(int a, int b)
{
	if(a > b) { return a; } else { return b; }
}

int min(int a, int b)
{
	if(a < b) { return a; } else { return b; }
}

int compar(const void* a, const void* b)
{
	return ( ((struct Phoneme*)a)->size_count - ((struct Phoneme*)b)->size_count );
}

void export_phones(void);
double** init_dtw_matrix(int signal_length, int phone_length, int w);
void* create_mfcc(void* argv);
void* create_clusters(void* argv);
void clean(void);
void gen_aoo(void);
void export_mfccs(void);
void export_clusters(void);
void average_mfccs(int num);
void export_for_jenks(char* p, int coef, float* data, int data_size);
void read_jenks(void);
void export_for_pca(void);
void normal_all(void);
void std_test_output(void);
void handle_argv(int argc, char* argv[]);
void cluster(void);
void pca(void);
void threaded_mfccs(void);
void min_max_values(void);
void mean_size_mfccs(struct Phoneme* phone);
void frame_variance(void);
void export_device(void);
void non_threaded_mfccs(void);

int seco(long time)
{
	return time - hour(time)*3600 - minu(time)*60;
}
int minu(long time)
{
	return (time - hour(time)*3600)/60;
}
int hour(long time)
{
	return time/3600;
}

int main(int argc, char* argv[])
{
	
	start = time(NULL);
	handle_argv(argc, argv);
	if((MALE + FEMALE + SPKR1 + SPKR1_NOSIL) > 1) {
		printf("Only one dataset may be chosen.\n If you wish to train both please do not input either...\n");
		exit(-1);
	}
	
	feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);
	printf("PAA_OP %d || PAA :: %d || BANKS :: %d || WINDOWS :: %d || FL :: %d || LIMIT :: %d || INTERVAL :: %d || NFFT : %d || MFCCS : %d || D :: %d || DD :: %d || COEFFS :: %.0f || EXTRA :: %d || CLUST :: %d || SVM :: %d || K :: %d || GROUP_K :: %d || VOICE_K :: %d || ITER :: %d\nZC :: %d || STE :: %d || ENTR :: %d || LARG_ENTR :: %d || LARG_STE :: %d\n", glbl_paa_op, glbl_paa, glbl_banks, glbl_window_width, glbl_frame_limit, glbl_dtw_window, glbl_interval_div, glbl_nfft, glbl_mfcc_num, DELTA, DELTA_DELTA, floor(glbl_banks * glbl_test_trunc), EXTRA, glbl_clust_num, SVM, glbl_k, glbl_group_k, glbl_voice_k, glbl_test_iter, glbl_zc_incr, glbl_ste_incr, glbl_entr_incr, glbl_larg_entr_incr, glbl_larg_ste_incr);
	
	printf(":: INITIALISING PHONEMES  ::  %02d:%02d:%02d  ::\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)));
	dtw_init();

	printf("::   TRAINING PHONEMES    ::  %02d:%02d:%02d  ::\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)));
	train();
	printf("::         TRAINED        ::  %02d:%02d:%02d  ::  %05d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), trained);

	printf("::     CREATING MFCCS     ::  %02d:%02d:%02d  ::\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)));
	if(THREAD) {
		threaded_mfccs();
	} else {
		non_threaded_mfccs();
	}

	
	struct sigaction action;
	action.sa_handler = interrupt_handler;
	sigemptyset(&action.sa_mask);
	action.sa_flags = SA_RESTART;
	if (sigaction(SIGINT, &action, NULL) == -1)
		printf("SIGINT handler failed\n");

	printf("::      CREATED MFCCS     ::  %02d:%02d:%02d  ::  %05d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), mfcced);

	printf("::    CREATING CLUSTERS   ::  %02d:%02d:%02d  ::  %05d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), clustered);

	// frame_variance();
	
	if(CLUST) {
		cluster();
	}

	if(PCA) {
		pca();
	}

	if(NORM) {
		printf("::      NORMALISING       ::  %02d:%02d:%02d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)));
		normal_all();
		printf("::    NORMALISING DONE    ::  %02d:%02d:%02d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)));
	}
	
	printf("::    CREATED CLUSTERS    ::  %02d:%02d:%02d  ::  %05d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), clustered);


	printf("::    TESTING PHONEMES    ::  %02d:%02d:%02d  ::\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)));
	for(int p = 0; p < glbl_test_iter; p++) {
		printf("Test iteration %d...\n", (p + 1));
		if(p == (glbl_test_iter - 1)) {
			// BOUNDS = 1;
		}
		test();
		float bytes = data_size();
		printf("::          BYTES         ::  %02d:%02d:%02d  ::  %.2fMB\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), bytes); 

		for(int i = 1; i < num_ph; i++) {
			for(int j = 0; j < phones[i]->size_count; j++) {
				if((phones[i]->correct[j] + phones[i]->error[j]) > 0.0 && phones[i]->size[j] != 0) {
					float per = phones[i]->correct[j] / (phones[i]->correct[j] + phones[i]->error[j]);
					if(per < 0.10F && phones[i]->reduced_count > 5 && phones[i]->size != 0) {
						phones[i]->size[j] = 0;
						phones[i]->reduced_count--;
						free(phones[i]->mfcc[j]);
					}
				}
			}
		}
		float new = data_size();
		printf("::          BYTES         ::  %02d:%02d:%02d  ::  %.2fMB\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), new); 
		if(!BOUNDS) {
			std_test_output();
		}
	}
	/* test(); */
	min_max_values();
	
	long double avg_time = total_test_time / total_dtw_tests;
	printf("::         TESTED         ::  %02d:%02d:%02d  ::  %05d:%Lfs\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), tested, avg_time);
	if(!BOUNDS) {
		std_test_output();
	}


	
	printf("::         LARGEST        ::  %02d:%02d:%02d  ::  %d:%f:%f:%d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), largest_mfcc, largest_value, smallest_value, largest_index);
	printf(":: DONE - MATRIX EXPORTED ::  %02d:%02d:%02d  ::\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start))); 
	printf("::       CLEANING UP      ::  %02d:%02d:%02d  ::\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)));
	
	if(CLUST) {
		export_clusters();
	}
	if(EXPORT) {
		export_device();
	}
	printf("::      EXPORTED MFCCS    ::  %02d:%02d:%02d  ::  %d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), exported);
	clean();
	
	
	return 0;
}

/** 
 * \brief Creates clusters for each phoneme
 * 
 * @param argv Not used
 *  
 * @return Not used
 *
 * This function is mostly unused due to the use of @file Jenks.py
 */
void* create_clusters(void* argv)
{
	pthread_mutex_lock(&lock);
	mask_sig();
	pthread_mutex_unlock(&lock);
	/* In testing, find cluster with   */
	int i = *(int*)argv;
	int total_frames = 0;		
	int trunc = floor(glbl_banks * glbl_test_trunc);
	
	
	int all_data_size = 0;
	for(int j = 0; j < phones[i]->size_count; j++) {
		for(int m = 0; m < phones[i]->size[j]; m++) {
			all_data_size++;
		}
	}
	
	if(all_data_size == 0) {
		perror("Can't have cluster with no data");
		exit(0);
	}
	
	float* all_data = (float*)calloc(all_data_size, sizeof(float));
	int t = 0;
	
	for(int j = 0; j < phones[i]->size_count; j++) {
		for(int m = 0; m < phones[i]->size[j]; m++) {
			// printf("i = %d | j = %d | m = %d | size = %d | count = %d\n", i, j, m, phones[i]->size[j], phones[i]->size_count);
			all_data[t] = phones[i]->mfcc[j][m];
			t++;
		}
	}	
	total_frames = t / trunc;
	float** data = (float**)malloc(trunc *  sizeof(float*));
	int n = 0;
	for(int p = 0; p < trunc; p++) {
		data[p] = (float*)calloc(total_frames, sizeof(float));
		for(int m = p; m < all_data_size; m+=trunc) {
			data[p][n] = all_data[m];
			n++;
		}
		n = 0;
	}
	phones[i]->clust = (struct Cluster**)malloc(sizeof(struct Cluster*) * trunc);
	time_t start = time(NULL);
	int final_sum = 0;
	pthread_mutex_lock(&lock);
	printf("::       CLUSTERING       ::  %-4s (%d) ->\n", phones[i]->index->name, all_data_size);
	pthread_mutex_unlock(&lock);
	float avggfv = 0;
	for(int p = 0; p < trunc; p++) {
		
		/* phones[i]->clust[p] = clust(data[p], total_frames, glbl_clust_num, 0, i); */
		/* final_sum += phones[i]->clust[p]->count; */
		/* avggfv += phones[i]->clust[p]->gfv; */
		
		pthread_mutex_lock(&lock);
		export_for_jenks(phones[i]->index->name, p, data[p], total_frames);
		pthread_mutex_unlock(&lock);
		free(data[p]);
	}
	// avggfv /= trunc;
	pthread_mutex_lock(&lock);
	clusters_left--;
	time_t end = time(NULL);
	printf(":: %-4s -> (%d) :: This took : %ld and achieved (%.2f) :: Clusters left : %d\n", phones[i]->index->name, final_sum, (end - start), avggfv, clusters_left);
	pthread_mutex_unlock(&lock);
		
	free(data);
	free(all_data);
	
	pthread_exit(NULL);
	return NULL;
}

/** 
 * \brief Creates all MFCCs for a given phoneme
 * 
 * @param argv The phoneme index to create MFCCs for
 * 
 * @return Not used
 *
 * Creates all MFCCs from all raw data signals for the phoneme index provided.
 */
void* create_mfcc(void* argv)
{
	pthread_mutex_lock(&lock);
	mask_sig();
	pthread_mutex_unlock(&lock);
	int i = *(int*)argv;
	// for(int i = 1; i < num_ph; i++) {
	// MFCCs
	phones[i]->mfcc = (float**)malloc(sizeof(float*) * phones[i]->size_count);
	phones[i]->mfcc_delta = (float**)malloc(sizeof(float*) * phones[i]->size_count);
	phones[i]->mfcc_delta_delta = (float**)malloc(sizeof(float*) * phones[i]->size_count);
	phones[i]->norm_mfcc = (float**)malloc(sizeof(float*) * phones[i]->size_count);
	// Features
	phones[i]->feats = (struct Feature_Set**)malloc(sizeof(struct Feature_Set*) * phones[i]->size_count);
	for(int p = 0; p < phones[i]->size_count; p++) {
		phones[i]->feats[p] = (struct Feature_Set*)malloc(sizeof(struct Feature_Set));
	}
	// Accuracy
	phones[i]->used = (int*)calloc(phones[i]->size_count, sizeof(int));
	phones[i]->correct = (float*)calloc(phones[i]->size_count, sizeof(int));
	phones[i]->error = (float*)calloc(phones[i]->size_count, sizeof(int));
	phones[i]->reduced_count = phones[i]->size_count;
	// Raw time
	if(RAW) {
		phones[i]->raw_sizes = (int*)calloc(phones[i]->size_count, sizeof(int));
	
		int uniq = 0;
		for(int m = 0; m < phones[i]->size_count; m++) {
			for(int n = 0; n < (uniq + 1); n++) {
				if(phones[i]->size[m] == phones[i]->raw_sizes[n]) {
					goto end;
				}
			}
			phones[i]->raw_sizes[uniq] = phones[i]->size[m];
			uniq++;
		end: ;
		}
		phones[i]->raw_count = uniq;
		phones[i]->raw_time = (float**)malloc(uniq * sizeof(float*));
		for(int m = 0; m < phones[i]->raw_count; m++) {
			for(int n = 0; n < phones[i]->size_count; n++) {
				if(phones[i]->raw_sizes[m] == phones[i]->size[n]) {
					phones[i]->raw_time[m] = (float*)calloc(phones[i]->size[n], sizeof(float));
					for(int p = 0; p < phones[i]->size[n]; p++) {
						phones[i]->raw_time[m][p] = (float)phones[i]->sequence[n][p];
					}
				}
			}
		}
	}
	for(int j = 0; j < phones[i]->size_count; j++) {
		float* signal = (float*)calloc(phones[i]->size[j], sizeof(float));
		short* short_signal = (short*)calloc(phones[i]->size[j], sizeof(short));
		float* zc_signal = (float*)calloc(phones[i]->size[j], sizeof(float));
		for(int n = 0; n < phones[i]->size[j]; n++) {
			signal[n] = phones[i]->sequence[j][n] / phones[i]->amounts[j];
			short_signal[n] = signal[n];
			zc_signal[n] = signal[n];
		}

		phones[i]->mfcc[j] = mfcc(signal, phones[i]->size[j], glbl_window_width, glbl_banks, glbl_paa_op);
		
		int new_size = mfcc_size(phones[i]->size[j]);
		int extra_size = floor(phones[i]->size[j] / glbl_window_width);
		phones[i]->feats[j]->coeffs = extra_size;
		if(EXTRA || (STE) || (ZC)) {
			phones[i]->feats[j]->zc = (float*)calloc(extra_size, sizeof(float));
			phones[i]->feats[j]->ste = (float*)calloc(extra_size, sizeof(float));
			phones[i]->feats[j]->kurtosis = (float*)calloc(extra_size, sizeof(float));
			phones[i]->feats[j]->entropy = (float*)calloc(extra_size, sizeof(float));
			for(int m = 0; m < extra_size; m++) {
				phones[i]->feats[j]->zc[m] = f_cross_rate(&zc_signal[m * glbl_window_width], glbl_window_width);
				phones[i]->feats[j]->ste[m] = short_time_energy(&short_signal[m * glbl_window_width], glbl_window_width, glbl_window_width);
				phones[i]->feats[j]->kurtosis[m] = kurtosis(&zc_signal[m * glbl_window_width], glbl_window_width);
				phones[i]->feats[j]->entropy[m] = log_entropy(&zc_signal[m * glbl_window_width], glbl_window_width);
			}
		}

		if(DELTA) {
			phones[i]->mfcc_delta[j] = delta(phones[i]->mfcc[j], new_size);
			if(DELTA_DELTA) {
				phones[i]->mfcc_delta_delta[j] = delta_delta(phones[i]->mfcc_delta[j], new_size);
			}
		}
			
		// memcpy(phones[i]->norm_mfcc[j], phones[i]->mfcc[j], new_size * sizeof(float));
							      
		phones[i]->size[j] = new_size;
		free(phones[i]->sequence[j]);
		free(short_signal);
		free(zc_signal);
		pthread_mutex_lock(&lock);
		mfcced++;
		if(mfcced % 1000 == 0)
			printf("::         MFCC P#        ::  %02d:%02d:%02d  ::  %05d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), mfcced);
		pthread_mutex_unlock(&lock);

				
	}
	free(phones[i]->sequence);
	// }

	if(AVG) {
		average_mfccs(i);
	} else if(MEAN_SIZE) {
		mean_size_mfccs(phones[i]);
	}

	

	pthread_exit(NULL);
	return NULL;
}

/** 
 * \brief Averages the MFCC coefficients for each phoneme.
 * 
 * @param num The phoneme index.
 */
void average_mfccs(int num)
{
	for(int i = 0; i < phones[num]->size_count; i++) {
		long double* sums = (long double*)calloc(phones[num]->size[i], sizeof(long double));
		
		int count = 0;
		if(phones[num]->size[i] == 0)
			continue;
		
		for(int j = 0; j < phones[num]->size_count; j++) {
			if( (phones[num]->size[j] == 0) || (phones[num]->size[j] != phones[num]->size[i]))
				continue;
			for(int m = 0; m < phones[num]->size[j]; m++) {
				sums[m] += phones[num]->mfcc[j][m];	
			}
			count++;
			if(i != j) {
				free(phones[num]->mfcc[j]);
				phones[num]->size[j] = 0;
			}
		}
		for(int f = 0; f < phones[num]->size[i]; f++) {
				phones[num]->mfcc[i][f] = sums[f] / count; 
			
		}

		free(sums);	
	}

	return;
}

/** 
 * \brief Initliases all phonemes found in @var p_codes
 *
 * All phonemes are initialised with the index, string code, MFCC array, raw data array, size array and MFCC count.
 */
void dtw_init(void)
{
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;
	ph_zc_max = calloc(num_ph, sizeof(float));
	ph_zc_min = calloc(num_ph, sizeof(float));
	ste_min = calloc(num_ph, sizeof(float));
	ste_max = calloc(num_ph, sizeof(float));
	for(int i = 1; i < num_ph; i++) {
		ph_zc_min[i] = INT_MAX;
	}

	phones = (struct Phoneme**)malloc(  (j + 1)  * sizeof(struct Phoneme*));
	num_ph = j;
	for(int i = 1; i < j; i++) {
		phones[i] = (struct Phoneme*)malloc(sizeof(struct Phoneme));
		phones[i]->size = (int*)malloc(sizeof(int));
		phones[i]->size[0] = 0;
		// Amounts array for averaging at the end
		phones[i]->amounts = (int*)malloc(sizeof(int));
		phones[i]->amounts[0] = 0;
		phones[i]->trained = 0;
		// Info about the phoneme
		phones[i]->index = (struct Ph_index*)malloc(sizeof(struct Ph_index));
		phones[i]->index->i = i; 
		phones[i]->index->name = strdup(p_codes[i]);
		phones[i]->size_count = 0;		
		if      (i >= 1 && i <= 6)   {
			memcpy(phones[i]->index->group, p_group[0],strlen(p_group[1]) + 1);
			phones[i]->index->group_i = 0;
			phones[i]->index->voice = 0;
		} else if (i == 7 || i == 8)   {
			memcpy(phones[i]->index->group, p_group[1], strlen(p_group[1]) + 1);
			phones[i]->index->group_i = 1;
			phones[i]->index->voice = 0;
		} else if (i >= 9 && i <= 15)  {
			memcpy(phones[i]->index->group, p_group[2], strlen(p_group[1]) + 1);
			phones[i]->index->group_i = 2;
			phones[i]->index->voice = 0;
		} else if (i >= 16 && i <= 18) {
			memcpy(phones[i]->index->group, p_group[3], strlen(p_group[1]) + 1);
			phones[i]->index->group_i = 3;
			phones[i]->index->voice = 1;
		} else if (i >= 19 && i <= 23) {
			memcpy(phones[i]->index->group, p_group[4], strlen(p_group[1]) + 1);
			phones[i]->index->group_i = 4;
			phones[i]->index->voice = 1;
		} else if (i >= 24 && i <= 37) {
			memcpy(phones[i]->index->group, p_group[5], strlen(p_group[1]) + 1);
			phones[i]->index->group_i = 5;
			phones[i]->index->voice = 1;
		} else {
			memcpy(phones[i]->index->group, p_group[6], strlen(p_group[1]) + 1);
			phones[i]->index->group_i = 6;
			phones[i]->index->voice = 2;
		}
	}
	 
	return;
} 

/** 
 * \brief A K-means function for use with @file jenks.py
 * 
 * @param signal The signal to be classified
 * @param signal_length The length of @param signal
 * @param phoneme The phoneme to be compared to
 * @param limit The window limit for DTW
 */
void dtw_clust(float** signal, int signal_length, struct Phoneme* phoneme, short limith)
{
	time_t started = time(NULL);
	phoneme->score = 0;
	double temp_last_min = 0,  last_min = 0;
	int phone_length = 0, w = 0; //  pos = 1;
	int p = 0;
	int mfcc_length = mfcc_size(signal_length);

	int trunc = floor((glbl_banks) * glbl_test_trunc);
	
	signal_length = mfcc_length / trunc;
	phone_length = signal_length;
	double score = 0;
	double cost = 0;
	
	int largest = max(signal_length, phone_length);
	float win = (float)glbl_dtw_window / 1000;
	w = floor(win * (float)largest);
	if(signal_length == 0) {
		printf("Exiting, no signal length() found\n");
		DTW_ERROR = NO_SIG_LEN;
		failed = 1;
		return;
	}
	if(phone_length == 0) {
		printf("Exiting, no sequence length found : %d : %d\n", phoneme->size[p], phone_length);
		failed = 1;
		return;
	}
	double** dtw_matrix = init_dtw_matrix(signal_length, signal_length, w);
	for(int i = 1; i <= signal_length - 1; i++) {
		for(int j = max(1, i-w); j <= min(phone_length - 1, i+w); j++) {
			for(int m = 0; m < trunc; m++) {
				double smallest_diff = DBL_MAX;
				for(int o = 0; o < phoneme->clust[m]->count; o++) {
					if(fabs(signal[0][(i * trunc) +  m] - phoneme->clust[m]->centroids[o]) < smallest_diff) {
						smallest_diff = fabs(signal[0][(i * trunc) +  m] - phoneme->clust[m]->centroids[o]);
						if(smallest_diff == 0)
							goto stop;
					}
				}
			stop:
				cost += pow(smallest_diff, 2);
			}
			temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
			last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
			dtw_matrix[i][j] = cost + last_min;
			score = dtw_matrix[i][j];
			cost = 0;
			// if(fabs(total_score + score) > best_so_far) { goto end; }
		}
	}

	phoneme->score = score;
		
	for(int i = 0; i < signal_length; i++) {
		free(dtw_matrix[i]);
	}
		
	free(dtw_matrix);
	time_t ended = time(NULL);
	total_test_time += (ended - started);
	total_dtw_tests++;
	return;
}

/** 
 * Similar to \fn dtw_clust() except it returns the result of DTW.
 */
double dtw_clust_result(float* signal, int signal_length, struct Phoneme* phoneme, int t, int clust, short limit)
{
	phoneme->score = 0;
	double temp_last_min = 0,  last_min = 0;
	int phone_length = 0, w = 0;
	int p = 0;
	int mfcc_length = mfcc_size(signal_length);
	double final_score = 0;
	int trunc = floor((glbl_banks) * glbl_test_trunc);
	
	signal_length = mfcc_length / trunc;
	phone_length = signal_length;
	double score = 0;
	double cost = 0;
	
	int largest = max(signal_length, phone_length);
	float win = (float)glbl_dtw_window / 1000;
	w = floor(win * (float)largest);
	if(signal_length == 0) {
		printf("Exiting, no signal length() found\n");
		DTW_ERROR = NO_SIG_LEN;
		failed = 1;
		return - 1;
	}
	if(phone_length == 0) {
		printf("Exiting, no sequence length found : %d : %d\n", phoneme->size[p], phone_length);
		failed = 1;
		return - 1;
	}
	double** dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
	for(int i = 1; i <= signal_length - 1; i++) {
		for(int j = max(1, i-w); j <= min(phone_length - 1, i+w); j++) {
			for(int m = 0; m < trunc; m++) {
				double smallest_diff = DBL_MAX;
				for(int o = 0; o < phoneme->clust[m]->count; o++) {
					if(fabs(signal[(i * trunc) +  m] - phoneme->clust[m]->centroids[o]) < smallest_diff) {
						smallest_diff = fabs(signal[(i * trunc) +  m] - phoneme->clust[m]->centroids[o]);
						if(smallest_diff == 0)
							goto stop;
						
					}
				}
			stop:
				cost += pow(smallest_diff, 2);
			}
			temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
			last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
			dtw_matrix[i][j] = cost + last_min;
			score = dtw_matrix[i][j];
			cost = 0;
			// if(fabs(total_score + score) > best_so_far) { goto end; }
		}
	}

	final_score = score;
		
	for(int i = 0; i < signal_length; i++) {
		free(dtw_matrix[i]);
	}
		
	free(dtw_matrix);
	return final_score;
}

/** 
 * \brief Performs DTW frame-by-frame for MFCCs
 * 
 * @param signal The input signal MFCC to be compared
 * @param signal_length The length of @param signal in frames
 * @param phoneme The phoneme to compare to 
 * @param limit The DTW window limit
 */
void dtw_frame(float** signal, int signal_length, struct Phoneme* phoneme, short limit)
{
	time_t started = time(NULL);
	phoneme->score = 0;
	double temp_last_min = 0,  last_min = 0;
	int phone_length = 0;
	int w = 0; //  pos = 1;
	int p = 0, smallest_diff = INT_MAX;
 	// float coef_diff = FLT_MAX;
	int mfcc_length = mfcc_size(signal_length);
	int to_try[phoneme->size_count];
	int loc = 0;
	for(int i = 0; i < phoneme->size_count; i++) {
		if(phoneme->size[i] == 0) {
			continue;
		}
		if((abs(phoneme->size[i] - mfcc_length) <= smallest_diff)) { // 
			smallest_diff = abs(phoneme->size[i] - mfcc_length);
			// coef_diff = fabs(signal[0][0] - phoneme->mfcc[i][0]);
			p = i;
		}
		if(smallest_diff == 0) {
			to_try[loc] = i;
			loc++;
		}
		
	}
	if(loc == 0) {
		to_try[0] = p;
		loc = 1;
	}

	double smallest = DBL_MAX;
	for(int m = 0; m < loc; m++) {
		p = to_try[m];
		double total_score = 0;
		int trunc = floor((glbl_banks) * glbl_test_trunc);	
		signal_length = mfcc_length / trunc;
		phone_length = phoneme->size[p] / trunc;
		double score = 0;
		double cost = 0;
		phoneme->used[p]++;
		/* An MFCC is amount * truncation - trunc = floor((incr / 2) * glbl_test_trunc)*/
		// w = max(limit, abs(signal_length-phone_length));
		int largest = max(signal_length, phone_length);
		float win = (float)glbl_dtw_window / 1000;
		w = floor(win * (float)largest);
		if(signal_length == 0) {
			printf("Exiting, no signal length() found\n");
			DTW_ERROR = NO_SIG_LEN;
			failed = 1;
			return;
		}
		if(phone_length == 0) {
			printf("Exiting, no sequence length found : %d : %d\n", phoneme->size[p], phone_length);
			failed = 1;
			return;
		}
		double** dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
		for(int i = 1; i <= signal_length - 1; i++) {
			for(int j = max(1, i-w); j <= min(phone_length - 1, i+w); j++) {
				//double** dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
				//for(int i = 1; i < signal_length - 1; i++) {
				//for(int j = max(1, i-w); j < min(phone_length - 1, i+w); j++) {
				for(int m = 0; m < trunc; m++) {
					cost += pow(fabs(signal[0][(i * trunc) +  m] - phoneme->mfcc[p][(j * trunc) +  m]), 2);
				}
				temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
				last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
				dtw_matrix[i][j] = cost + last_min;
				score = dtw_matrix[i][j];
				cost = 0;
				// if(fabs(total_score + score) > best_so_far) { goto end; }
			}
		}
		total_score += score;
		
		for(int i = 0; i < signal_length; i++) {
			free(dtw_matrix[i]);
		}
		
		free(dtw_matrix);

		if(DELTA) {
			float* signal_delta = delta(signal[0], mfcc_length);
			if(NORM) {
				normalise_delta(signal_delta, mfcc_length);
			}
			dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
			score = 0;
			for(int i = 1; i < signal_length - 1; i++) {
				for(int j = max(1, i-w); j < min(phone_length - 1, i+w); j++) {
					for(int m = 0; m < trunc; m++) {
						cost += fabs(signal_delta[(i * trunc) +  m] - phoneme->mfcc_delta[p][(j * trunc) +  m]);
					}

					temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
					last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
					dtw_matrix[i][j] = cost + last_min;
					score = dtw_matrix[i][j];
					cost = 0;
					// if(fabs(score) > best_so_far) { goto end; }
				}
			}

			total_score += score;
		
			for(int i = 0; i < signal_length; i++) {
				free(dtw_matrix[i]);
			}
		
			free(dtw_matrix);
			free(signal_delta); 
		}
		if(DELTA_DELTA) {
			float* signal_delta_delta = delta(signal[0], mfcc_length);
			if(NORM) {
				normalise_delta_delta(signal_delta_delta, mfcc_length);
			}
			dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
			score = 0;
			for(int i = 1; i < signal_length - 1; i++) {
				for(int j = max(1, i-w); j < min(phone_length - 1, i+w); j++) {
					for(int m = 0; m < trunc; m++) {
						cost += fabs(signal_delta_delta[(i * trunc) +  m] - phoneme->mfcc_delta_delta[p][(j * trunc) +  m]);
					}
					temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
					last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
					dtw_matrix[i][j] = cost + last_min;
					score = dtw_matrix[i][j];
					cost = 0;
					// if(fabs(score) > best_so_far) { goto end; }
				}
			}


			total_score += score;
// end:

			for(int i = 0; i < signal_length; i++) {
				free(dtw_matrix[i]);
			}
			free(dtw_matrix); 
			free(signal_delta_delta);
		}
		if(fabs(total_score) < smallest || m == 0) {
			smallest = fabs(total_score);
			phoneme->score = smallest;
		}
	}

	time_t ended = time(NULL);
	total_test_time += (ended - started);
	total_dtw_tests++;
	return;
}

/**
 * \brief Similar to \fn dtw_frame() except the result is returned.
 */
double dtw_frame_result(float* signal, int signal_length, struct Phoneme* phoneme, int p, short limit)
{
	
	phoneme->score = 0;
	double temp_last_min = 0,  last_min = 0;
	int phone_length = 0;
	int w = 0;
	int mfcc_length = mfcc_size(signal_length);
	double final_score = 0;
	int trunc = floor((glbl_banks) * glbl_test_trunc);
	signal_length = mfcc_length / trunc;
	phone_length = phoneme->size[p] / trunc;
	
	double score = 0;
	double cost = 0;
	phoneme->used[p]++;
	int largest = max(signal_length, phone_length);
	float win = (float)glbl_dtw_window / 1000;
	w = floor(win * (float)largest);
	// printf("Window :: %d :: (%d * %f)\n", w, largest, win);
	if(signal_length == 0) {
		printf("Exiting, no signal length() found\n");
		DTW_ERROR = NO_SIG_LEN;
		failed = 1;
		return -1;
	}
	if(phone_length == 0) {
		printf("Exiting, no sequence length found : %d : %d\n", phoneme->size[p], phone_length);
		failed = 1;
		return -1;
	}
	int diff = 1, start = 1;

	double** dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
	for(int i = start; i <= signal_length - 1; i++) {
		for(int j = max(start, i-w); j <= min(phone_length - diff, i+w); j++) {
			for(int m = 0; m < trunc; m++) {
				cost += fabs(signal[(i * trunc) +  m] - phoneme->mfcc[p][(j * trunc) +  m]); // * weight);
			}
			temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
			last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
			dtw_matrix[i][j] = cost + last_min;
			score = dtw_matrix[i][j];
			cost = 0;
		}
	}
	
	final_score += score;	
	for(int i = 0; i < signal_length; i++) {
		free(dtw_matrix[i]);
	}
	free(dtw_matrix);

	if(DELTA) {
		float* signal_delta = delta(signal, mfcc_length);
		if(NORM) {
			normalise_delta(signal_delta, mfcc_length);
		}
		dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
		score = 0;
		for(int i = start; i <= signal_length - diff; i++) {
			for(int j = max(start, i-w); j <= min(phone_length - diff, i+w); j++) {
				for(int m = 0; m < trunc; m++) {
					cost += fabs(signal_delta[(i * trunc) +  m] - phoneme->mfcc_delta[p][(j * trunc) +  m]);
				}

				temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
				last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
				dtw_matrix[i][j] = cost + last_min;
				score = dtw_matrix[i][j];
				cost = 0;
				// if(fabs(score) > best_so_far) { goto end; }
			}
		}

		final_score += score;
		
		for(int i = 0; i < signal_length; i++) {
			free(dtw_matrix[i]);
		}
		
		free(dtw_matrix);
		free(signal_delta); 
	}
	if(DELTA_DELTA) {
		float* signal_delta_delta = delta(signal, mfcc_length);
		if(NORM) {
			normalise_delta_delta(signal_delta_delta, mfcc_length);
		}
		dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
		score = 0;
		for(int i = start; i <= signal_length - diff; i++) {
			for(int j = max(start, i-w); j <= min(phone_length - diff, i+w); j++) {
				for(int m = 0; m < trunc; m++) {
					cost += fabs(signal_delta_delta[(i * trunc) +  m] - phoneme->mfcc_delta_delta[p][(j * trunc) +  m]);
				}
				temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
				last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
				dtw_matrix[i][j] = cost + last_min;
				score = dtw_matrix[i][j];
				cost = 0;
				// if(fabs(score) > best_so_far) { goto end; }
			}
		}


		final_score += score;
// end:

		for(int i = 0; i < signal_length; i++) {
			free(dtw_matrix[i]);
		}
		free(dtw_matrix); 
		free(signal_delta_delta);
	}

	return final_score;
}

/**
 * \brief Similar to \fn dtw_frame_result() except for group classication
 *
 * This method can classify a phoneme's group using MFCCs, zero cross, short time energy, or deltas, depending on the set parameter.
 */
double dtw_frame_result_group(float** signal, int signal_length, struct Phoneme* phoneme, int p, short limit)
{
	
	phoneme->score = 0;
	double temp_last_min = 0,  last_min = 0, final_score = 0;
	int phone_length = 0;
	int w = 0;
	int mfcc_length = mfcc_size(signal_length);
	int trunc = floor((glbl_banks) * glbl_test_trunc);
	int s_length = signal_length;
	int signal_coeffs = floor(s_length / glbl_window_width);
	
	signal_length = mfcc_length / trunc;
	phone_length = phoneme->size[p] / trunc;
	double score = 0;
	double cost = 0;
	phoneme->used[p]++;
	int largest = max(signal_length, phone_length);
	float win = (float)glbl_dtw_window / 1000;
	double** dtw_matrix;
	w = floor(win * (float)largest);
	if(signal_length == 0) {
		printf("Exiting, no signal length() found\n");
		DTW_ERROR = NO_SIG_LEN;
		failed = 1;
		return -1;
	}
	if(phone_length == 0) {
		printf("Exiting, no sequence length found : %d : %d\n", phoneme->size[p], phone_length);
		failed = 1;
		return -1;
	}

	if(DELTA) {
		float* signal_delta = delta(signal[0], mfcc_length);
		if(NORM) {
			normalise_delta(signal_delta, mfcc_length);
		}
		dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
		score = 0;
		for(int i = start; i <= signal_length - 1; i++) {
			for(int j = max(start, i-w); j <= min(phone_length - 1, i+w); j++) {			       
				for(int m = 0; m < trunc; m++) {
					cost += fabs(signal_delta[(i * trunc) +  m] - phoneme->mfcc_delta[p][(j * trunc) +  m]);
				}

				temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
				last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
				dtw_matrix[i][j] = cost + last_min;
				score = dtw_matrix[i][j];
				cost = 0;
				// if(fabs(score) > best_so_far) { goto end; }
			}
		}

		final_score += score;
		
		for(int i = 0; i < signal_length; i++) {
			free(dtw_matrix[i]);
		}
		
		free(dtw_matrix);
		free(signal_delta);
	}
	if(DELTA_DELTA) {
		float* signal_delta_delta = delta(signal[0], mfcc_length);
		if(NORM) {
			normalise_delta_delta(signal_delta_delta, mfcc_length);
		}
		dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
		score = 0;
		for(int i = start; i <= signal_length - 1; i++) {
			for(int j = max(start, i-w); j <= min(phone_length - 1, i+w); j++) {
				for(int m = 0; m < trunc; m++) {
					cost += fabs(signal_delta_delta[(i * trunc) +  m] - phoneme->mfcc_delta_delta[p][(j * trunc) +  m]);
				}
				temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
				last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
				dtw_matrix[i][j] = cost + last_min;
				score = dtw_matrix[i][j];
				cost = 0;
				// if(fabs(score) > best_so_far) { goto end; }
			}
		}


		final_score += score;
// end:

		for(int i = 0; i < signal_length; i++) {
			free(dtw_matrix[i]);
		}
		free(dtw_matrix); 
		free(signal_delta_delta);
	}
	
	if(GROUP && (ZC || STE) && !(DELTA || DELTA_DELTA)) {
		w = max(limit, abs(signal_coeffs-phoneme->feats[p]->coeffs));
		dtw_matrix = init_dtw_matrix(signal_coeffs, phoneme->feats[p]->coeffs, w);
		score = 0;
		for(int i = 1; i < signal_coeffs - 1; i++) {
			for(int j = max(1, i-w); j < min(phoneme->feats[p]->coeffs - 1, i+w); j++) {
				if(ZC) {
					cost += fabs(phoneme->feats[p]->zc[j] - signal[1][i]);
				} else if(STE){
					cost += fabs(phoneme->feats[p]->ste[j] - signal[2][i]);
				} else if(KURT) {
					cost += fabs(phoneme->feats[p]->kurtosis[j] - signal[3][i]);
				} else if(ENTR) {
					cost += fabs(phoneme->feats[p]->entropy[j] - signal[4][i]);
				}
				temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
				last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
				dtw_matrix[i][j] = cost + last_min;
				score = dtw_matrix[i][j];
				cost = 0;
			}
		}

		final_score += score;
		
		for(int i = 0; i < signal_coeffs; i++) {
			free(dtw_matrix[i]);
		}
		
		free(dtw_matrix);
		return final_score;
	}

	
	
	return final_score;
}

/**
 * \brief Similar to \fn dtw_frame_result_group() except it uses raw time signals
 *
 * This method can classify a phoneme's group using raw time signals.
 */
double dtw_frame_result_group_time(float* signal, int signal_length, struct Phoneme* phoneme, int p, short limit)
{
	
	phoneme->score = 0;
	double temp_last_min = 0,  last_min = 0, final_score = 0;
	int phone_length = 0;
	int w = 0;	
	signal_length = signal_length;
	phone_length = phoneme->raw_sizes[p];
	double score = 0;
	double cost = 0;
	// phoneme->used[p]++;
	int largest = max(signal_length, phone_length);
	float win = (float)glbl_dtw_window / 1000;
	double** dtw_matrix;
	w = floor(win * (float)largest);
	if(signal_length == 0) {
		printf("Exiting, no signal length() found\n");
		DTW_ERROR = NO_SIG_LEN;
		failed = 1;
		return -1;
	}
	if(phone_length == 0) {
		printf("Exiting, no sequence length found : %d : %d : %d\n", phoneme->raw_sizes[p], p, phoneme->raw_count);
		failed = 1;
		return -1;
	}

	if(DELTA) {
		dtw_matrix = init_dtw_matrix(signal_length, phone_length, w);
		score = 0;
		for(int i = start; i <= signal_length - 1; i++) {
			for(int j = max(start, i-w); j <= min(phone_length - 1, i+w); j++) {
			       
				cost += fabs(signal[i] - phoneme->raw_time[p][j]);
				temp_last_min = fmin(dtw_matrix[i-1][j], dtw_matrix[i][j-1]);
				last_min = fmin(temp_last_min, dtw_matrix[i-1][j-1]);
				dtw_matrix[i][j] = cost + last_min;
				score = dtw_matrix[i][j];
				cost = 0;
				// if(fabs(score) > best_so_far) { goto end; }
			}
		}

		final_score += score;
		
		for(int i = 0; i < signal_length; i++) {
			free(dtw_matrix[i]);
		}
		
		free(dtw_matrix);
	}
	
	
	return final_score;
}

/** 
 * \brief Initalises a zeroed matrix for use in DTW
 * 
 * @param signal_length The length of the test's sequence MFCC in frames
 * @param phone_length The length of the phoneme's sequence MFCC in frames
 * @param w The windowing limit
 * 
 * @return The zero initialised matrix
 */
double** init_dtw_matrix(int signal_length, int phone_length, int w)
{
	double** dtw_matrix = (double**)malloc(signal_length *  sizeof(double*));
	if(dtw_matrix == NULL)
	{
		fprintf(stderr, "out of memory\n");
		exit(-1);
	}
	for(int i = 0; i < signal_length; i++) {
		dtw_matrix[i] = (double*)malloc(phone_length * sizeof(double));
		if(dtw_matrix[i] == NULL)
			printf("Failed to malloc 'dtw_matrix[i]'\n");
	}
	for(int i = 0; i < signal_length; i++) {
		for(int j = 0; j < phone_length; j++) {
			dtw_matrix[i][j] = DBL_MAX;
		}
	}
	dtw_matrix[0][0] = 0;

	for(int i = 1; i <= signal_length - 1; i++) {
		for(int j = max(1, i-w); j <= min(phone_length - 1, i+w); j++) {
			//for(int i = 1; i < signal_length; i++) {
			//for(int j = max(1, i-w); j < min(phone_length, i+w); j++) {
			dtw_matrix[i][j] = 0;
		}
	}

	return dtw_matrix;
}

/** 
 * \brief Returns the length of an array terminated with SHRT_MAX. 
 * 
 * @param array The array to find the length of
 * 
 * @return The length of @param array
 *
 * No longer used in any function as the size is passed instead of finding it.
 */
size_t length(short* array)
{
	int size = 0, i = 0;
	while(array[i] != SHRT_MAX) {
		size++;
		i++;
	}
	return size;
}

int is_sil_mean(short* array, int length)
{
	double total = 0, mean = 0; //, std = 0;
	for(int i = 0; i < length - 2; i++) {
		total += array[i];
	}
	mean = total / (length - 2);

	if(mean > max_sil_mean || mean < min_sil_mean) {
		return 1;
	}
	
	float ste = short_time_energy(array, length - 2, glbl_window_width);
	float* signal = (float*)malloc(sizeof(float) * (length - 2));
	for(int i = 0; i < length - 2; i++) {
		signal[i] = array[i];
	}
	int zc = stavg_cross_rate_no_overlap(signal, length - 2);
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;
	for(int i = 1; i < j; i++) {
		// Get the appropriate ste and zc to compare
		if(strcmp(phones[i]->index->name, "h#") == 0 ||
		   strcmp(phones[i]->index->name, "epi") == 0 ||
		   strcmp(phones[i]->index->name, "pau") == 0) {
			// compare the ste and zc
			if(!(zc < ph_zc_max[i] && zc > ph_zc_min[i] &&
			     ste < ste_max[i] && ste > ste_min[i])) {
				// if these pass then compare against the sil min-max values
				return 1;
				
			}
		}
	}

	free(signal);
	return 0;
}

int is_sil(short* array, int length)
{
	double sum = 0, total = 0, mean = 0; //, std = 0;
	for(int i = 0; i < length - 2; i++) {
		total += array[i];
	}
	mean = total / (length - 2);
	for(int i = 0; i < length - 2; i++) {
		sum += pow( (array[i] - mean), 2);
	}
	// std = sqrt(sum / (length - 2));
	float ste = short_time_energy(array, length - 2, glbl_window_width);
	float* signal = (float*)malloc(sizeof(float) * (length - 2));
	for(int i = 0; i < length - 2; i++) {
		signal[i] = array[i];
	}
	int zc = stavg_cross_rate_no_overlap(signal, length - 2);
	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;
	for(int i = 1; i < j; i++) {
		// Get the appropriate ste and zc to compare
		if(strcmp(phones[i]->index->name, "h#") == 0 ||
		   strcmp(phones[i]->index->name, "epi") == 0 ||
		   strcmp(phones[i]->index->name, "pau") == 0) {
			// compare the ste and zc
			if(!(zc < ph_zc_max[i] && zc > ph_zc_min[i] &&
			     ste < ste_max[i] && ste > ste_min[i])) {
				// if these pass then compare against the sil min-max values
				return 1;
				
			}
		}
	}
	for(int i = 0; i < length - 2; i++) {
		//int Z_score = (array[i] - mean) / std;
		//if(Z_score < 2 && Z_score > -2) {
		//	continue;
		//}
		if(array[i] > max_sil || array[i] < min_sil) {
			free(signal);
			return 1;
		}
	}
	free(signal);
	return 0;
}

int is_sil_ste_chunk(short* array, int length)
{

	float amount = floor(length / glbl_window_width), ste = 0, count = 0;
	float* signal = (float*)calloc(length, sizeof(float));
	for(int i = 0; i < length; i++) {
		signal[i] = array[i];
	}
	for(int i = 0; i < (length - glbl_window_width) ; i+=glbl_window_width) {
		ste = f_short_time_energy(&signal[i], glbl_window_width, glbl_window_width);
		if(ste > 2500)
			count++;
		if( (count / amount) > 0.1)
			return 1;
	}
	return 0;

}

int is_sil_ste(short* array, int length)
{
	float* signal = (float*)calloc(length, sizeof(float));
	for(int i = 0; i < length; i++) {
		signal[i] = array[i];
	}
	float ste = f_short_time_energy(signal, length, glbl_window_width);
	if(ste > 30000) {
		free(signal);
		return 1;
	}
	if(ste > 10000) {
		free(signal);
		return 2;
	}
		
	free(signal);
	return 0;
}

int is_sil_zc(short* array, int length)
{
	float n = 0;
	float zc = 0;
	for(int i = 0; i < length - glbl_window_width; i+=glbl_window_width) {
		zc += favg_cross_rate(&array[i], glbl_window_width);
		n++;
	}
	zc /= n;
	// printf("zc :: %f\n", zc);
	if(zc > 0.3) {
		return 1;
	}
	if(zc > 0.2) {
		return 2;
	}

	return 0;
}

int is_sil_flat(short* array, int length)
{
	float* sig = (float*)calloc(length, sizeof(float));
	for(int i = 0; i < length; i++) {
		sig[i] = array[i];
	}
	float flat = flatness(sig, length); 
	// printf("zc :: %f\n", zc);
	if(flat < min_sil_flt || flat > max_sil_flt) {
		return 1;
	}
	free(sig);
	return 0;
}

int is_sil_db(short* array, int length)
{

	for(int i = 0; i < length - 1; i++) {
		float amplitude = abs(array[i]) / 32767;
		if(amplitude == 0)
			amplitude = FLT_EPSILON;
		float dB = 20 * log10(amplitude);
		if(dB > 15 || dB < -15) {
			return 1;
		}
	}

	return 0;
}

/** 
 * \brief Cleans up any used memory before closing the program.
 * 
 */
void clean(void)
{

	int j = 1;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}

	num_ph = j;
	for(int i = 1; i < j; i++) {

		// Inital sequence
		for(int j = 0; j < phones[i]->size_count; j++) {
			if(phones[i]->size[j] == 0)
				continue;
			free(phones[i]->mfcc[j]);
			if(DELTA) { 
				free(phones[i]->mfcc_delta[j]);
			}
			if(DELTA_DELTA) {
				free(phones[i]->mfcc_delta_delta[j]);
			}
			if(EXTRA || ZC || STE) {
				free(phones[i]->feats[j]->zc);
				free(phones[i]->feats[j]->ste);
				free(phones[i]->feats[j]->kurtosis);
				free(phones[i]->feats[j]->entropy);
				free(phones[i]->feats[j]);
			}
			free(phones[i]->feats[j]);
			
			// free(phones[i]->norm_mfcc[j]);
			// free(phones[i]->aao_i);
		
		}
		if(RAW) {
			for(int j = 0; j < phones[i]->raw_count; j++) {
				free(phones[i]->raw_time[j]);
			}
			free(phones[i]->raw_time);
		}
		free(phones[i]->feats);
		free(phones[i]->correct);
		free(phones[i]->error);
		free(phones[i]->used);
		free(phones[i]->size);
		free(phones[i]->amounts);
		free(phones[i]->mfcc);
		free(phones[i]->mfcc_delta);
		free(phones[i]->mfcc_delta_delta);
		free(phones[i]->norm_mfcc);
		free(phones[i]->index->name);
		free(phones[i]->index);
		free(phones[i]);
	}
	free(phones);
	return;
}

/** 
 * \brief Closes all threads in use
 * 
 */
void cancel_threads(void)
{
	for(int i = 0; i < thread_count; i++) {
		if(threads[i] != 0) {
			pthread_kill(threads[i], SIGKILL);
		}
	}
}

/** 
 * \brief Handles a signal interrupt from the user
 * 
 * @param signo The signal flag
 *
 * When a CTRL+C (SIGINT) is sent to the program the used memory is cleaned and the current stats are displayed before exiting
 */
void interrupt_handler(int signo)
{
	if(signo == SIGINT) {
		perror("\nSIGNIT! - Cleaning Up!\n");
		if(!BOUNDS) {
			float per = 0;
			per = ((float)correct / ((float)group + (float)correct + (float)fails)) * 100.0;
			long double avg_time = total_test_time / total_dtw_tests;
			printf("Current correct :: %.2f%%\nTotal tested :: %d\nAverage test time :: %Lfs\n", per, tested, avg_time);
			printf("Elapsed time    ::  %02d:%02d:%02d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)));
		}
		if(GROUP) { 
			float g_corr = 0, g_err = 0, g_per = 0;
			for(int i = 0; i < 7; i++) {
				for(int j = 0; j < 7; j++) {
					if(j == i) {
						g_corr += group_matrix[i][j];
					} else {
						g_err += group_matrix[i][j];
					}
				}
		
			}
			g_per = (g_corr / (g_corr + g_err)) * 100;
			printf("::          GROUP         ::            ::  %.2f%% \n", g_per);
			float v_corr = 0, v_err = 0, v_per = 0, sv_corr = 0, sv_err = 0, sv_per = 0;
			for(int i = 0; i < 3; i++) {
				for(int j = 0; j < 3; j++) {
					if(j == i) {
						v_corr += voice_matrix[i][j];
						sv_corr += sil_v_matrix[i][j];
					} else {
						v_err += voice_matrix[i][j];
						sv_err += sil_v_matrix[i][j];
					}
				}
		
			}
			sv_err = 0;
			v_per = (v_corr / (v_corr + v_err)) * 100;
			printf("::          VOICE         ::            ::  %.2f%%:%.2f%% \n", v_per, sv_per);
		}
		if(BOUNDS) {
			long double end_wer = ((ins+subs+delts)/(refs));
			printf("::           WER          ::  %02d:%02d:%02d  ::  %.2Lf\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), end_wer);
		}
		cancel_threads();
		clean();
		
		exit(EXIT_FAILURE);
	}

	return;
}

/** 
 * \brief Masks any signals for additional threads other than main
 * 
 */
void mask_sig(void)
{
	sigset_t mask;
	sigemptyset(&mask); 
        sigaddset(&mask, SIGRTMIN+3); 
                
        pthread_sigmask(SIG_BLOCK, &mask, NULL);
        
}

/** 
 * \brief Exports all MFCCs in use to files in the main directory
 * 
 */
void export_mfccs(void)
{
	
	char* base = "../MFCCs/";
	char sub_fold[1024];
	struct dirent *next_file;
	char filepath[1024];
	for(int i = 1; i < num_ph; i++) {
		sprintf(sub_fold, "%s%s/", base, phones[i]->index->name);
		DIR* theFolder = opendir(sub_fold);
		if(theFolder == NULL) {
			printf("Couldn't open folder :: %s\n", sub_fold);
			continue;
		}
		while ( (next_file = readdir(theFolder)) != NULL )
		{
			if(strcmp(next_file->d_name, "../MFCCs/") == 0)
				continue;
			// build the path for each file in the folder
			sprintf(filepath, "%s", sub_fold);
			sprintf(filepath, "%s/", next_file->d_name);
			remove(filepath);
		}
		closedir(theFolder);
	}
	
	
	char filename[256];
	FILE* m_fp = fopen("../MFCCs/config.txt", "w");
	fprintf(m_fp, "banks %d\ntrunc %f\nwidth %d\noverlap %d\n", glbl_banks, glbl_test_trunc, glbl_window_width, glbl_interval_div);
	fclose(m_fp);
	
	int trunc = floor(glbl_banks * glbl_test_trunc);
	for(int i = 1; i < num_ph; i++) {
		for(int j = 0; j < phones[i]->size_count; j++) {
			if(phones[i]->size[j] != 0) {
				sprintf(filename, "%s/%s/", base, phones[i]->index->name);
				struct stat st = {0};
				if (stat(filename, &st) == -1) {
					mkdir(filename, 0700);
				}
				sprintf(filename, "%s/%s/%d_%f.txt", base, phones[i]->index->name, phones[i]->size[j], phones[i]->mfcc[j][0]);
				FILE* fp = fopen(filename, "wb");
				int n = 0;
				for(int m = 0; m < phones[i]->size[j]; m++) {
					fprintf(fp, "%f ", phones[i]->mfcc[j][m]);
					n++;
					if(n == (trunc)) {
						fprintf(fp, "\n");
						n = 0;
					}
				}
				exported++;
				fclose(fp);
			}
		}
	}

	return;
}

/** 
 * \brief Reads in the Jenks centroids generated by @file jenks.py
 * 
 */
void read_jenks(void)
{
	char* path = "../Jenks/";
	char filename[512];
	char* end_ptr;
	int trunc = floor(glbl_banks * glbl_test_trunc);
	for(int i = 1; i < num_ph; i++) {
		int n = 0;
		char line[256];
		char *p;
		phones[i]->clust = (struct Cluster**)malloc(sizeof(struct Cluster*) * trunc);
		for(int t = 0; t < trunc; t++) {
			phones[i]->clust[t] = (struct Cluster*)malloc(sizeof(struct Cluster));
			sprintf(filename, "%s%s/res/%d.jen", path, phones[i]->index->name, t);
			FILE* fp = fopen(filename, "r");
			while(fgets(line, sizeof(line), fp)) {
				n++;
			}
			// printf("file :: %s\n", filename);
			phones[i]->clust[t]->centroids = (float*)calloc(n, sizeof(float));
			n = 0;
			fseek(fp, SEEK_SET, 0);
			while(fgets(line, sizeof(line), fp)) {
				p = strtok(line,"\n");
				// printf("p :: %s\n", p);
			        phones[i]->clust[t]->centroids[n] = strtol(p, &end_ptr, 10);
				n++;
			}
			phones[i]->clust[t]->count = n;
		}
	}

	return;
}

/** 
 * \brief Exports all MFCCs to files so that @file jenks.py can process them
 * 
 * @param p The phoneme code to be exported
 * @param coef The coefficient number to be exported
 * @param data The MFCC data for @param p and @param coef
 * @param data_size The length of the @param data array
 */
void export_for_jenks(char* p, int coef, float* data, int data_size)
{
	
	char* base = "../Jenks/";
	char sub_fold[512];
	sprintf(sub_fold, "%s%s/", base, p);
	struct stat st1 = {0};
	if (stat(sub_fold, &st1) == -1) {
		mkdir(sub_fold, 0700);
	}
		
	char filename[256];
	
	sprintf(filename, "%s/%s/", base, p);
	struct stat st = {0};
	if (stat(filename, &st) == -1) {
		mkdir(filename, 0700);
	}
	sprintf(filename, "%s/%s/%d.txt", base, p, coef);
	FILE* fp = fopen(filename, "w");
	for(int m = 0; m < data_size; m++) {
		fprintf(fp, "%f,", data[m]);
	}
	exported++;
	fclose(fp);


	return;
}

/** 
 * \brief Exports all MFCCs in use for use with @file pca.py
 * 
 */
void export_for_pca(void)
{
	
	char* base = "../PCA_Data/";
	char sub_fold[512];
	sprintf(sub_fold, "%s_Pca/", base);
	struct stat st1 = {0};
	if (stat(sub_fold, &st1) == -1) {
		mkdir(sub_fold, 0700);
	}
		
	char filename[1024];
	int trunc = floor(glbl_banks * glbl_test_trunc);
	
	FILE* fp = NULL;
	sprintf(filename, "%s/data_.txt", sub_fold);
	fp = fopen(filename, "a");
       
	for(int i = 1; i < num_ph; i++) {
		for(int j = 0; j < phones[i]->size_count; j++) {
			for(int m = 0; m < phones[i]->size[j] / trunc; m++) {
				for(int t = 0; t < trunc; t++) {
					fprintf(fp, "%f,", phones[i]->mfcc[j][(m * trunc) + t]);			
				}
				fprintf(fp, "%s\n", phones[i]->index->name);
			}
		}
	}
	fclose(fp);
	exported++;
	return;
}

/** 
 * \brief Exports all clusters in use
 * 
 */
void export_clusters(void)
{
	
	char* base = "../Clusters/";
	char sub_fold[1024];
	struct dirent *next_file;
	char filepath[2048];
	for(int i = 1; i < num_ph; i++) {
		sprintf(sub_fold, "%s%s/", base, phones[i]->index->name);
		DIR* theFolder = opendir(sub_fold);
		if(theFolder == NULL) {
			printf("Couldn't open folder :: %s\n", sub_fold);
			continue;
		}
		while ( (next_file = readdir(theFolder)) != NULL )
		{
			if(strcmp(next_file->d_name, "../Clusters/") == 0)
				continue;
			// build the path for each file in the folder
			sprintf(filepath, "%s%s/", sub_fold, next_file->d_name);
			remove(filepath);
		}
		closedir(theFolder);
	}
	
	
	char filename[1024];
	FILE* m_fp = fopen("../Clusters/config.txt", "w");
	fprintf(m_fp, "banks %d\ntrunc %f\nwidth %d\noverlap %d\n", glbl_banks, glbl_test_trunc, glbl_window_width, glbl_interval_div);
	fclose(m_fp);
	
	int trunc = floor(glbl_banks * glbl_test_trunc);
	for(int i = 1; i < num_ph; i++) {
		for(int m = 0; m < trunc; m++) {
			sprintf(filename, "%s/%s/", base, phones[i]->index->name);
			struct stat st = {0};
			if (stat(filename, &st) == -1) {
				mkdir(filename, 0700);
			}
			sprintf(filename, "%s/%s/%d.txt", base, phones[i]->index->name, m);
			FILE* fp = fopen(filename, "wb");
			int n = 0;
			for(int j = 0; j < phones[i]->clust[m]->count; j++) {
				fprintf(fp, "%f ", phones[i]->clust[m]->centroids[j]);
				n++;
				if(n == 16) {
					fprintf(fp, "\n");
					n = 0;
				}
			}
			
			exported++;
			fclose(fp);
		}
		
	}

	return;
}

/** 
 * \brief Exports all MFCCs and config files for use with the physical device
 *
 * This method produces a folder to be copied onto the device's SD card, all required information will be loaded
 * by the device
 */
void export_device(void)
{
	
	char* base = "../device/";
	DIR* base_fold = opendir(base);
	if(base_fold == NULL) {
		printf("No /device/ folder, making one now...\n");
		mkdir(base, 0700);
	}
	closedir(base_fold);
	char sub_fold[256];
	struct dirent *next_file;
	char filepath[512];
	for(int i = 1; i < num_ph; i++) {
		sprintf(sub_fold, "%s%s/", base, phones[i]->index->name);
		DIR* theFolder = opendir(sub_fold);
		if(theFolder == NULL) {
			printf("Couldn't open folder :: %s - making one...\n", sub_fold);
			mkdir(sub_fold, 0700);
			closedir(theFolder);
			theFolder = opendir(sub_fold);
		}
		while ( (next_file = readdir(theFolder)) != NULL )
		{
			if(strcmp(next_file->d_name, "../base/") == 0)
				continue;
			// build the path for each file in the folder
			sprintf(filepath, "%s%s/", sub_fold, next_file->d_name);
			remove(filepath);
		}
		closedir(theFolder);
	}
	
	
	char filename[1024];
	FILE* m_fp = fopen("../device/config.conf", "w");
	fprintf(m_fp, "banks %d\ntrunc %f\nwidth %d\noverlap %d\nframe_limit %d", glbl_banks, glbl_test_trunc, glbl_window_width, glbl_interval_div,
		glbl_frame_limit);
	fclose(m_fp);
	
	// int trunc = floor(glbl_banks * glbl_test_trunc);
	for(int i = 1; i < num_ph; i++) {
		sprintf(filename, "%s/%s/", base, phones[i]->index->name);
		struct stat st = {0};
		if (stat(filename, &st) == -1) {
			mkdir(filename, 0700);
		}
		sprintf(filename, "%s/%s/%s.conf", base, phones[i]->index->name, phones[i]->index->name);
		FILE* fp = fopen(filename, "wb");
		fprintf(fp, "%d \n", phones[i]->reduced_count);
		for(int j = 0; j < phones[i]->size_count; j++) {
			if(phones[i]->size[j] == 0)
				continue;
			fprintf(fp, "%d ", phones[i]->size[j]);
		}
		fprintf(fp, " \n");
		int n = 0;
		int max = 0;
		for(int j = 0; j < phones[i]->size_count; j++) {
			if(phones[i]->size[j] > max) {
				max = phones[i]->size[j];
			}
		}
		int* amounts = calloc((max+1), sizeof(int));
		for(int j = 0; j < phones[i]->size_count; j++) {
			if(phones[i]->size[j] == 0) {
				continue;
			}
			sprintf(filename, "%s/%s/%d_%d.phn", base, phones[i]->index->name, phones[i]->size[j], amounts[phones[i]->size[j]]);
			
			amounts[phones[i]->size[j]]++;
			fp = freopen(filename, "wb", fp);
			for(int m = 0; m < phones[i]->size[j]; m++) {
				fprintf(fp, "%f ", phones[i]->mfcc[j][m]);
				n++;
			}
			fprintf(fp, "\n");
			exported++;
			sprintf(filename, "%s/%s/%s.conf", base, phones[i]->index->name, phones[i]->index->name);
			fp = freopen(filename, "ab", fp);
			fprintf(fp, "%d ", amounts[phones[i]->size[j]]);
		}
		sprintf(filename, "%s/%s/%s.conf", base, phones[i]->index->name, phones[i]->index->name);
		fp = freopen(filename, "ab", fp);
		fprintf(fp, " \n");
		fclose(fp);
		free(amounts);
	}

	return;
}

/** 
 * \brief Normalises all MFCCs and deltas
 * 
 */
void normal_all(void)
{

	for(int i = 1; i < num_ph; i++) {
		for(int j = 0; j < phones[i]->size_count; j++) {
			update_mfcc_norm(phones[i]->mfcc[j], phones[i]->size[j]);
		}
	}
	for(int i = 1; i < num_ph; i++) {
		for(int j = 0; j < phones[i]->size_count; j++) {
			normalise_mfcc(phones[i]->mfcc[j], phones[i]->size[j]);
		}
	}
	if(DELTA) {
		for(int i = 1; i < num_ph; i++) {
			for(int j = 0; j < phones[i]->size_count; j++) {
				update_deltas_norm(phones[i]->mfcc_delta[j], phones[i]->size[j]);
			}
		}
		for(int i = 1; i < num_ph; i++) {
			for(int j = 0; j < phones[i]->size_count; j++) {
				normalise_delta(phones[i]->mfcc_delta[j], phones[i]->size[j]);
			}
		}
	}
	if(DELTA_DELTA) {
		for(int i = 1; i < num_ph; i++) {
			for(int j = 0; j < phones[i]->size_count; j++) {
				update_delta_deltas_norm(phones[i]->mfcc_delta_delta[j], phones[i]->size[j]);
			}
		}
		for(int i = 1; i < num_ph; i++) {
			for(int j = 0; j < phones[i]->size_count; j++) {
				normalise_delta_delta(phones[i]->mfcc_delta_delta[j], phones[i]->size[j]);
			}
		}
	}

	return;
}

/** 
 * \brief Produces standard (non-boundary) tests output to the user
 * 
 */
void std_test_output(void)
{
	
	float per = 0;
	
	per = ((float)correct / ((float)group + (float)correct + (float)fails)) * 100.0;
	float sil_per = 0;
	sil_per = ((float)sil_corr / ((float)correct + (float)sil_corr)) * 100.0;
	float no_win = 0;
	if(no_window != 0) {
		no_win = no_window / tested;
	} else {
		no_win = 0;
	}
	float zero_per = 0;
	if((zero_fails + fails) >= 0) {
		zero_per = (float)zero_fails / ((float)zero_fails + (float)fails);
	}
	
	if(GROUP) { 
		float g_corr = 0, g_err = 0, g_per = 0;
		for(int i = 0; i < 7; i++) {
			for(int j = 0; j < 7; j++) {
				if(j == i) {
					g_corr += group_matrix[i][j];
				} else {
					g_err += group_matrix[i][j];
				}
			}
		
		}
		g_per = (g_corr / (g_corr + g_err)) * 100;
		printf("::          GROUP         ::            ::  %.2f%% \n", g_per);
		float v_corr = 0, v_err = 0, v_per = 0, sv_corr = 0, sv_err = 0, sv_per = 0;
		for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) {
				if(j == i) {
					v_corr += voice_matrix[i][j];
					sv_corr += sil_v_matrix[i][j];
				} else {
					v_err += voice_matrix[i][j];
					sv_err += sil_v_matrix[i][j];
				}
			}
		
		}
		sv_err = 0;
		if((v_corr + v_err) == 0) {
			v_per = 0;
		} else {
			v_per = (v_corr / (v_corr + v_err)) * 100;
		}
		printf("::          VOICE         ::            ::  %.2f%%:%.2f%% \n", v_per, sv_per);
	}
	
	printf("::         RESULTS        ::            ::  %.2f%% : %.2f%% : %.2f : %.2f%%\n", per, sil_per, no_win, zero_per);

	return;
}

/** 
 * \brief Handles all input arguments provided by the user from the command line
 * 
 * @param argc The number of arguments provided
 * @param argv The arguments as strings
 */
void handle_argv(int argc, char* argv[])
{


	if(argc > 1) {
		char* end_ptr;
		for(int i = 1; i < argc; i++) {
			if(strcmp(argv[i], "paa") == 0) {
				glbl_paa = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "window") == 0) {
				glbl_window_width = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "banks") == 0) {
				glbl_banks = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "paa_op") == 0) {
				glbl_paa_op = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "clust_num") == 0) {
				glbl_clust_num = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "dtw_window") == 0) {
				glbl_dtw_window = strtol(argv[i + 1], &end_ptr, 10); i++;
				//glbl_dtw_window /= 1000;				
			} else if(strcmp(argv[i], "interval_div") == 0) {
				glbl_interval_div = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "knn") == 0) {
				glbl_k = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "test_iter") == 0) {
				glbl_test_iter = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "nfft") == 0) {
				glbl_nfft = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "mfccs") == 0) {
				glbl_mfcc_num = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "zc_incr") == 0) {
				glbl_zc_incr = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "entr_incr") == 0) {
				glbl_entr_incr = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "ste_incr") == 0) {
				glbl_ste_incr = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "neg_incr") == 0) {
				glbl_neg_incr = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "larg_entr_incr") == 0) {
				glbl_larg_entr_incr = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "larg_ste_incr") == 0) {
				glbl_larg_ste_incr = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "larg_incr") == 0) {
				glbl_larg_incr = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "group_k") == 0) {
				glbl_group_k = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "voice_k") == 0) {
				glbl_voice_k = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "frame_limit") == 0) {
				glbl_frame_limit = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "DELTA") == 0) {
				DELTA = 1;
			} else if(strcmp(argv[i], "DELTA_DELTA") == 0) {
				DELTA_DELTA = 1;
			} else if(strcmp(argv[i], "GROUP") == 0) {
				GROUP = 1;
			} else if(strcmp(argv[i], "GRAM") == 0) {
				GRAM = 1;
			} else if(strcmp(argv[i], "NORM") == 0) {
				NORM = 1;
			} else if(strcmp(argv[i], "SIMPLE") == 0) {
				SIMPLE = 1;
			} else if(strcmp(argv[i], "MALE") == 0) {
				MALE = 1;
			} else if(strcmp(argv[i], "FEMALE") == 0) {
				FEMALE = 1;
			} else if(strcmp(argv[i], "RAW") == 0) {
				RAW = 1;
			} else if(strcmp(argv[i], "SPKR1") == 0) {
				SPKR1 = 1;
			} else if(strcmp(argv[i], "SPKR1_NOSIL") == 0) {
				SPKR1_NOSIL = 1;
			} else if(strcmp(argv[i], "MEAN_SIZE") == 0) {
				MEAN_SIZE = 1;
			} else if(strcmp(argv[i], "VOICED") == 0) {
				VOICED = 1;
			} else if(strcmp(argv[i], "THREAD") == 0) {
				THREAD = 1;
			} else if(strcmp(argv[i], "EXPORT") == 0) {
				EXPORT = 1;
			} else if(strcmp(argv[i], "trunc") == 0) {
				glbl_test_trunc = strtol(argv[i + 1], &end_ptr, 10); i++;
			} else if(strcmp(argv[i], "EXTRA") == 0) {
				EXTRA = 1;
			} else if(strcmp(argv[i], "KNN") == 0) {
				KNN = 1;
			} else if(strcmp(argv[i], "CLUST") == 0) {
				CLUST = 1;
			} else if(strcmp(argv[i], "BOUNDS") == 0) {
				BOUNDS = 1;
			} else if(strcmp(argv[i], "PCA") == 0) {
				PCA = 1;
			} else if(strcmp(argv[i], "LOG_E") == 0) {
				LOG_E = 1;
			} else if(strcmp(argv[i], "AVG") == 0) {
				AVG = 1;
			} else if(strcmp(argv[i], "STE") == 0) {
				STE = 1;
			} else if(strcmp(argv[i], "ZC") == 0) {
				ZC = 1;
			} else if(strcmp(argv[i], "STRT") == 0) {
				STRT = 1;
			} else if(strcmp(argv[i], "WHITE") == 0) {
				WHITE = 1;
			} else if(strcmp(argv[i], "CAR") == 0) {
				CAR = 1;
			} else if(strcmp(argv[i], "CAFE") == 0) {
				CAFE = 1;
			} else {
				printf("Unknown argument : %s\nExiting...\n", argv[i]);
				exit(-1);
			}
		}
	}
	glbl_test_trunc = (float)glbl_test_trunc / (float)glbl_banks;
	if(glbl_paa == 1) {
		glbl_paa_op = 1;
	}
	SIMPLE = 1;
	if(CAFE || CAR || STRT || WHITE) {
		printf("Noise type has been selected, selecting SPKR1 dataset...\n");
		SPKR1 = 1;
	}
	int data_set = (SPKR1 + SPKR1_NOSIL + FEMALE + MALE);
	if(data_set > 1) {
		printf("More than one dataset chosen, please only choose one\n");
		exit(0);
	}
	return;
}

/** 
 * \brief Produces clusters for MFCCs using @file jenks.py
 * 
 */
void cluster(void)
{
	setvbuf(stdout, NULL, _IONBF, 0);
	pthread_t* threads_clust = calloc(num_ph, sizeof(pthread_t));
	int* values_clust = (int*)calloc(num_ph, sizeof(int));
	for(int i = 1; i < num_ph; i++)
		clusters_left++;
	
	for(int i = 1; i < num_ph; i++) {
		values_clust[i] = i;
		int err = pthread_create(&threads_clust[i], NULL, create_clusters, (void *)&values_clust[i]);
		if(err) {
			printf("Cluster thread error, exiting...\n");
			exit(-1);
		}
	}
	for(int j = 1; j < num_ph; j++) {
		int err = pthread_join(threads_clust[j], NULL);
		if (err) {
			printf("Thread join failed\n");
			exit(-1);
		}
	}
	free(threads_clust);
	free(values_clust);
	setvbuf(stdout, NULL, _IOLBF, 0);
	printf("Calling upon a serpent...\n");
	char run[512];
	sprintf(run, "python3 jenks.py %d", glbl_clust_num);
	system(run);
	printf("Reading Jenks...\n");
	read_jenks();
	printf("Done!\n");

	return;
}

/** 
 * \brief Exports the MFCCs for use with @file pca.py then cleans and exits.
 * 
 */
void pca(void)
{

	export_for_pca();
	// int trunc = floor(glbl_banks * glbl_test_trunc);
	printf("Calling upon a serpent...\n");
	char run[512];
	// sprintf(run, "python3 pca.py %d", trunc);
	system(run);
	printf("Done, got look for yourself...\n");
	clean();
	exit(0);

	return;
}

/** 
 * \brief Generates a thread for each phoneme to produce all MFCCs for that phoneme
 * 
 */
void threaded_mfccs(void)
{
	pthread_mutex_init(&lock, NULL);
	pthread_t* threads = calloc(num_ph, sizeof(pthread_t));
	int* values = (int*)calloc(num_ph, sizeof(int));

	for(int i = 1; i < num_ph; i++) {
		values[i] = i;
		int err = pthread_create(&threads[i], NULL, create_mfcc, (void *)&values[i]);
		if(err) {
			printf("Cluster thread error, exiting...\n");
			exit(-1);
		}
	}
	for(int j = 1; j < num_ph; j++) {
		int err = pthread_join(threads[j], NULL);
		if (err) {
			printf("Thread join failed\n");
			exit(-1);
		}
	}
	free(threads);
	free(values);

	return;
}

/** 
 * \brief Generates a thread for each phoneme to produce all MFCCs for that phoneme, only one thread is used.
 * 
 */
void non_threaded_mfccs(void)
{
	pthread_mutex_init(&lock, NULL);
	pthread_t* threads = calloc(num_ph, sizeof(pthread_t));
	int* values = (int*)calloc(num_ph, sizeof(int));

	for(int i = 1; i < num_ph; i++) {
		values[i] = i;
		int err = pthread_create(&threads[i], NULL, create_mfcc, (void *)&values[i]);
		if(err) {
			printf("Cluster thread error, exiting...\n");
			exit(-1);
		}
		err = pthread_join(threads[i], NULL);
		if (err) {
			printf("Thread join failed\n");
			exit(-1);
		}
	}
	free(threads);
	free(values);

	return;
}

/** 
 * \brief Assigns the minimum and maximum values found in all MFCCs
 * 
 */
void min_max_values(void)
{

	for(int i = 1; i < num_ph; i++) {		
		for(int j = 0; j < phones[i]->size_count; j++) {
			if(phones[i]->size[j] == 0)
				continue;
			if(phones[i]->size[j] > largest_mfcc) {
				largest_mfcc = phones[i]->size[j];
				largest_index = i;
			}
			for(int m = 0; m < phones[i]->size[j]; m++) {
				if(phones[i]->mfcc[j][m] > largest_value) {
					largest_value = phones[i]->mfcc[j][m];
				}
				if(phones[i]->mfcc[j][m] < smallest_value) {
					smallest_value = phones[i]->mfcc[j][m];
				}
			}
		}
	}

	return;
}

/** 
 * \brief Averages all MFCCs of the same size by phoneme
 * 
 * @param phone The phoneme to average MFCCs for
 */
void mean_size_mfccs(struct Phoneme* phone)
{
	int amount = 0;
	for(int i = 0; i < phone->size_count; i++) {
		if(phone->size[i] <= 0)
			continue;
		for(int j = 0; j < phone->size_count; j++) {
			if(j == i)
				continue;
			if(phone->size[i] == phone->size[j]) {
				for(int m = 0; m < phone->size[j]; m++) {
					phone->mfcc[i][m] += phone->mfcc[j][m];
				}
				amount++;
				phone->size[j] = 0;
			}
			
		}
		if(amount == 0)
			continue;
		for(int m = 0; m < phone->size[i]; m++) {
			phone->mfcc[i][m] /= amount;
		}
		amount = 0;	
	}
	
	return;
}

/** 
 * \brief Calculates and outputs the variance among frames of all MFCCs.
 * 
 */
void frame_variance(void)
{
	int lng = 0;
	int trunc = floor(glbl_banks * glbl_test_trunc);
	for(int i = 1; i < num_ph; i++) {
		for(int j = 0; j < phones[i]->size_count; j++) {
			if(phones[i]->size[j] == 0)
				continue;
			if(phones[i]->size[j] / trunc > lng) {
				lng = phones[i]->size[j] / trunc;
			}
		}
	}
	for(int p = 0; p < lng; p++) {
		long double sum = 0, variance = 0, mean = 0, n = 0;
		for(int i = 1; i < num_ph; i++) {
			for(int j = 0; j < phones[i]->size_count; j++) {
				if(phones[i]->size[j] / trunc > p) {
					long double temp_sum = 0;
					for(int t = 0; t < trunc; t++) {
						temp_sum = phones[i]->mfcc[j][(p * trunc) + t];
						
					}
					sum += (temp_sum / trunc);
					n++;
				}
			}
		}
		if(n == 0)
			continue;
		mean = (sum / n);
		for(int i = 1; i < num_ph; i++) {
			for(int j = 0; j < phones[i]->size_count; j++) {
				if(phones[i]->size[j] / trunc > p) {
					for(int t = 0; t < trunc; t++) {
						variance += pow((phones[i]->mfcc[j][(p * trunc) + t] - mean), 2);
					}
				}
			}
		}
		variance /= n;
		printf("Variance for frame %d is %Lf\n", p, variance);
	}

	return;
}




