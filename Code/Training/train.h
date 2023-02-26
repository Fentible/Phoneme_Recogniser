#ifndef TRAIN_H
#define TRAIN_H

#include <dirent.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <sys/types.h>

#include "../Clustering/cluster.h"
#include "../Dynamic_Time_Warping/dtw.h"
#include "../Misc/realloc.h"
#include "../Testing/test.h"
#include "../Feature_Extraction/paa.h"
#include "../Feature_Extraction/mfcc.h"
#include "../Seperation/cross_rate.h"
#include "../Seperation/ste.h"

#define TRAIN 0
#define TEST  1

struct wav_file {
	short* seq;
	int length;
	int offset;
};

double cubic_interpolate(short y0, short y1, short y2, short y3, double mu);
void train(void);
struct wav_file* read_wav(FILE* fp);
void allocate_ph(FILE* fp, short* wav, unsigned char t_t);
short* train_ph(int new, short* sequence, struct Phoneme* ph); 
short* init_new_phone(struct Phoneme* phone, short* sequence, int new);
short* resize(short* shorter, size_t s, size_t l);
int mfcc_size(int signal_length);
int frame_amount(int signal_length);
long double* ld_resize(long double* shorter, size_t s, size_t l);
void update_sil_zc(short* array, int length, char* filename);
void update_sil_ste(short* array, int length, char* filename);
float flatness(float* chunk, int length);

extern struct Phoneme** phones;
extern char* p_codes[];
extern char* p_group[];
extern int num_ph;
extern int prev_ph;

extern int limit_changed;

extern float* ph_zc_max;
extern float* ph_zc_min;
extern long double* zc_st_avg_no_oc;
extern float* ste_min;
extern float* ste_max;

extern int max_sil;
extern int min_sil;
extern float max_sil_dB;
extern float min_sil_dB;
extern float max_sil_zc;
extern float min_sil_zc;
extern float max_sil_ste;
extern float min_sil_ste;
extern float max_sil_flt;
extern float min_sil_flt;
extern float max_sil_mean;
extern float min_sil_mean;

#endif
