#ifndef TEST_H
#define TEST_H

#include <math.h>
#include "../Clustering/cluster.h"
#include "../Dynamic_Time_Warping/dtw.h"
#include "../Training/train.h"
#include "../Seperation/cross_rate.h"
#include "../Seperation/bounds.h"
#include "../Clustering/knn.h"

void export_results(char* ph_code);
void export_results_aao(char* ph_code);
void test(void);
void test_phoneme(short* h, int signal_length, char* p);
void test_phoneme_aao(short* h, long start_end[2], char* p);
void test_phoneme_pca(short* h, long start_end[2], char* p);
float reduce_data(void);
float data_size(void);
void means(void);

extern int correct;
extern int fails;
extern int group;
extern int sil_corr;
extern int zero_fails;

extern float* per_correct;

extern int* shorted;
extern int* removed;
extern int shted;
extern int rmved;

extern long double refs;
extern long double ins;
extern long double delts;
extern long double subs;

extern int group_matrix[7][7];
extern int voice_matrix[3][3];
extern int sil_v_matrix[3][3];

#endif

