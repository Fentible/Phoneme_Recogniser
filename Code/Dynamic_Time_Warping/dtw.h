#ifndef DTW_H
#define DTW_H

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <fenv.h>
#include <errno.h>
#include <pthread.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ftw.h>

#include "../Clustering/cluster.h"

struct Ph_index {
	int i;
	char group[5];
	int group_i;
	char* name;
	int voice;
};

struct Feature_Set {
	float* entropy;
	float* kurtosis;
	float* zc;
	float* ste;
	int coeffs;
	// float flatness;
	// float skewness;
	// float spread;	
}; 

struct Phoneme {
	struct Ph_index* index;
	struct Cluster** clust;
	struct Feature_Set** feats;
	float** mfcc_delta;
	float** mfcc_delta_delta;
	float** norm_mfcc;
	double score;
	float** mfcc;
	long double** sequence;
	// int* coeffs;
	int* size;
	int* used;
	int* amounts;
	int size_count;
	int reduced_count;
	long double trained;
	float* correct;
	float* error;
	float** raw_time;
	int* raw_sizes;
	int raw_count;
} ph;

#include "../Training/train.h"
#include "../Misc/realloc.h"
#include "../Testing/test.h"
#include "../Feature_Extraction/paa.h"
#include "../Feature_Extraction/mfcc.h"
#include "../Feature_Extraction/delta.h"

#define PTHREAD_CANCELED ((void *) -1)

// Global options - may be best to make them local
extern int glbl_banks;
extern int glbl_window_width;
extern int glbl_paa_op;
extern int glbl_paa;
extern int glbl_dtw_window;
extern int glbl_interval_div;
extern int glbl_nfft;
extern int glbl_mfcc_num;
extern int no_window;
extern int glbl_clust_num;
extern int glbl_k;
extern int glbl_test_iter;
extern int glbl_group_k;
extern int glbl_voice_k;
extern int glbl_frame_limit;

extern int glbl_zc_incr;
extern int glbl_ste_incr;
extern int glbl_entr_incr;
extern int glbl_neg_incr;
extern int glbl_larg_incr;
extern int glbl_larg_entr_incr;
extern int glbl_larg_ste_incr;

extern int DELTA;
extern int DELTA_DELTA;
extern float glbl_test_trunc;

#define NO_SIG_LEN 1
extern int DTW_ERROR;
extern int SIMPLE;
extern int MALE;
extern int FEMALE;
extern int EXTRA;
extern int CLUST;
extern int BOUNDS;
extern int Z_ZC;
extern int ONE;
extern int LOG_E;
extern int KNN;
extern int GROUP;
extern int TO_RUN;
extern int GRAM;
extern int NORM;
extern int ZC;
extern int STE;
extern int VOICED;
extern int ENTR;
extern int KURT;
extern int SPKR1;
extern int SPKR1_NOSIL;
extern int WHITE;
extern int STRT;
extern int CAR;
extern int CAFE;
extern int SPLIT_DATA;
extern int SVM;

// Testing externs 
extern float* aaomfcc;
extern int aaomfcc_length;
extern double against_all[];

// Output globals
extern double best_so_far;
extern short failed;
extern int trained;
extern int tested;
extern time_t start;

extern long double total_test_time;
extern long double total_dtw_tests;

extern pthread_t* threads;
extern int thread_count;

void dtw(float* signal, int signal_length, struct Phoneme* phoneme, short limit);
void dtw_init();
size_t length(short* array);
void export_phone(struct Phoneme* phoneme);
void dtw_test(float* signal, float* sequence, int signal_length, int signal_length2, short limit);
int seco(long time);
int minu(long time);
int hour(long time);
int is_sil(short* array, int length);
void dtw_frame(float** signal, int signal_length, struct Phoneme* phoneme, short limit);
int min(int a, int b);
int max(int a, int b);
void interrupt_handler (int signo);
int is_sil_db(short* array, int length);
int is_sil_mean(short* array, int length);
void dtw_clust(float** signal, int signal_length, struct Phoneme* phoneme, short limit);
void mask_sig(void);
double dtw_frame_result(float* signal, int signal_length, struct Phoneme* phoneme, int p, short limit);
double dtw_clust_result(float* signal, int signal_length, struct Phoneme* phoneme, int t, int clust, short limit);
double dtw_frame_result_group(float** signal, int signal_length, struct Phoneme* phoneme, int p, short limit);
double dtw_frame_result_group_time(float* signal, int signal_length, struct Phoneme* phoneme, int p, short limit);
int is_sil_zc(short* array, int length);
int is_sil_ste(short* array, int length);
int is_sil_flat(short* array, int length);
int is_sil_ste_chunk(short* array, int length);


#endif
