/**
 * \file cluster.c
 * \brief A naive Jenks clustering class - \file jenks.py is now used instead
 */
#include "cluster.h"

/**
 * \fn breaks
 * \brief finds and return the breaks within a data set based on the change in value
 * @return an array of breaks in the dataset 
 * @data the data to be processed
 * @length the length of the data array
 * @amount the amount of breaks to be returned
 * @offset if set to 1 a random offset is set to the iteration of the generated breaks 
 */
float* breaks(float* data, int length, int amount, int offset);

/** 
 * \fn mean
 * \brief computes and returns the mean of a dataset - used to compute the GVF
 * @data the array of values
 * @length length of the given array
 */
float mean(float* data, int length);

/**
 * \fn remove_dups
 * \brief moves any duplicates to the end of the array and returns the size of the array up to the duplicate values
 * @data the array to be processed
 * @length the length of the given array
 *
 * No data is actually removed in the process of this function, only moved
 */
int remove_dups(float* data, int length);

/** 
 * \fn remove_zeros
 * \brief moves any zero values to the end of the array and returns size of the array up to the zero values
 * @data the array to be processed
 * @length the length of the given array
 * 
 * No data is actually removed in the process of this function, only moved
 */ 
int remove_zeros(float* data, int length);

int fltcomp(const void * a, const void * b);

/**
 * \fn threadproc
 * A simple function each thread calls to output its current progress every 5 minutes
 */
void *threadproc(void *arg);

/**
 * \fn clust
 * \brief The main function which takes in an array of values and returns a cluster structure
 * @data the array of values to be processed
 * @length the length of the given array
 * @amount the lower bound of clusters to be returned
 * @offset whether to apply a random offset to the breaks function
 * @num the phoneme index
 * 
 * In this current state offset is set based on the progress of the function and amount
 * has an upper bound set to amount ^ 2. 
 * 
 * \warning This function is random and can be slow. Please use the jenks.py option instead.
 */
struct Cluster* clust(float* data, int length, int amount, int offset, int num);

/** 
 * \fn SDCM
 * \brief calculates the squared deviations from the class means used for calculating the GVF
 * @data the array to be processed
 * @length the length of the given array
 * @return the SDCM value of the array 
 */
float SDCM(float* data, int length);

int tries[64] = {0};         /* The number of tries that a thread has completed before the \var EX[num] limit has been reset */
float EX[64] = {0};          /* The current GVF limit for each thread to exit */
float tries_mult[64] = {1};  /* How many times a GVF limit has been reset for each thread */

static pthread_mutex_t lock = PTHREAD_MUTEX_INITIALIZER;

int fltcomp(const void * a, const void * b) {
	float fa = *(const float*) a;
	float fb = *(const float*) b;
	if (fa < fb) return -1;
	if (fa > fb) return 1;
	return 0;
}

void *threadproc(void *arg)
{
	int num = *(int*)arg;
	while(1) 
	{
		sleep(5 * 60);
		printf("... Still going ... [%f - (%f)] - [%f] - (%d)\n", tries[num] * tries_mult[num], tries_mult[num], EX[num], num);
	}
	pthread_exit(NULL);
}

struct Cluster* clust(float* data, int length, int amount, int offset, int num)
{
	EX[num] = 0.8;
	tries[num] = 0;
	int min = amount;
	amount = amount * amount;
	pthread_mutex_lock(&lock);
	pthread_t tid;
	int error = pthread_create(&tid, NULL, &threadproc, &num);
	if(threads == NULL) {
		threads = (pthread_t*)malloc(sizeof(pthread_t));
		threads[thread_count] = tid;
		thread_count++;
	} else {
		threads = realloc(threads, sizeof(pthread_t) * (thread_count + 1));
		threads[thread_count] = tid;
		thread_count++;
	}
	pthread_mutex_unlock(&lock);
	if(error) {
		printf("Sorry, error making cluster timer\n");
	}
	while(amount >= length) {
		amount /= 4;
	}
	struct Cluster* result = (struct Cluster*)malloc(sizeof(struct Cluster));
	length = length - 1;
	length = remove_dups(data, length);
	qsort(data, length, sizeof(float), fltcomp);
	int k = amount; 
	float GVF = 0;
	while(1) {
		int size = k; 
		float* brks = breaks(data, length - 1, k, offset);
		qsort(brks, k - 1, sizeof(float), fltcomp);
		result->centroids = (float*)calloc(k, sizeof(float));
		for(int i = 0; i < (k - 1); i++) {
			result->centroids[i] = (brks[i] + brks[i + 1]) / 2;
		}
		result->count = remove_dups(result->centroids, size);
		result->count = remove_zeros(result->centroids, result->count);
		result->values = (float**)malloc(result->count * sizeof(float*));
		for(int i = 0; i < result->count; i++) {
			result->values[i] = (float*)calloc(1, sizeof(float));
		}
		result->sizes = (int*)calloc(result->count, sizeof(int));
		for(int i = 0; i < result->count; i++) {
			result->sizes[i] = 0;
		}
		float diff = FLT_MAX;
		int pos = 0;
		for(int i = 0; i < length - 1; i++) {
			for(int j = 0; j < result->count; j++) {
				if(fabs(result->centroids[j] - data[i]) < diff) {
					diff = fabs(result->centroids[j] - data[i]);
					pos = j;
				}
			}
			result->values[pos][result->sizes[pos]] = data[i];
			result->sizes[pos]++;
			result->values[pos] = f_realloc(result->values[pos], result->sizes[pos] * sizeof(float));
			result->values[pos][result->sizes[pos]] = 0;
			diff = FLT_MAX;
		}
		float SDAM = 0, SDC = 0;
		for(int i = 0; i < result->count; i++) {
			SDC += SDCM(result->values[i], result->sizes[i]);
		}
		SDAM = SDCM(data, length); 
		GVF = (SDAM - SDC) / SDAM;
		offset = 1;
		tries[num]++;
		if(GVF >= EX[num]) {
			qsort(result->centroids, result->count, sizeof(float), fltcomp);
			free(brks);
			pthread_mutex_lock(&lock);
			for(int i = 0; i < thread_count; i++) {
				if(tid == threads[i]) {
					threads[i] = 0;
				}
			}
			pthread_cancel(tid);
			pthread_mutex_unlock(&lock);
			for(int i = 0; i < result->count; i++) {
				free(result->values[i]);
			}
			free(result->sizes);
			free(result->values);
			result->gfv = EX[num];
			return result;
		} else {
			// If the GVF is not above the threshold
			free(result->centroids);
			for(int i = 0; i < result->count; i++) {
				free(result->values[i]);
			}
			free(result->sizes);
			free(result->values);
			// Reset the GVF threshold and reduce the amount of clusters to be returned
			if(k <= min) {
				k = (amount - 1);
				EX[num] = 1 - ( log(tries[num]) / (log(100) / log(2)) );
				tries_mult[num]++;
			} else {
				k--;
			}
			// If below the minimum desired GVF reset to a high GVF
			if(EX[num] <= 0.25) {
				EX[num] = 0.8;
			}
			// If the thread has tried too many times and no good GVF can be found then accept any GVF
			if(tries_mult[num] >= 1000) {
				EX[num] = 0.0;
			}
		}
		free(brks);
	}
	
	
	return result;
}

float* breaks(float* data, int length, int amount, int offset)
{
	float* break_s = (float*)calloc(length, sizeof(float));
	float* result = (float*)calloc(amount, sizeof(float));
	// Iterate of the data and produce an array of value change percentages
	for(int i = 0; i < (length - 1); i++) {
		if(fabs(data[i]) == 0) {
			break_s[i] = 0;
		}
		else {
			break_s[i] = 100 * (data[i + 1] - data[i]) / fabs(data[i]); 
		}
	}

	// Get a random value for the offset - unix timestamp is too slow and random/urandom is too slow and blocks
	if(offset == 1) {
		struct timespec ts;
		if (timespec_get(&ts, TIME_UTC) == 0) {
			printf("Error getting seed, returning\n");
			return NULL;
		}
		srandom(ts.tv_nsec ^ ts.tv_sec);
		offset = random() % (length - 1); 
	}
	float largest = 0;
	int place = 0;
	// Get the breaks with the highest percentage
	for(int j = 0; j < (amount - 1); j++) {
		for(int i = offset; i < (length - 1); i++) {
			if(break_s[i] > largest) {
				largest = break_s[i];
				place = i;
			}
		}
		result[j] = data[place]; 
		break_s[place] = 0;
		largest = 0;
	}
	free(break_s);
	return result;

}

float mean(float* data, int length)
{
	if(length == 0)
		return 0;
	int sum = 0;
	float result = 0;
	for(int i = 0; i < length; i++) {
		sum += data[i];
	}
	result = sum / length;
	return result;

}
	
float SDCM(float* data, int length)
{
	if(length == 0)
		return 0;
	float diff = 0, sum = 0;
	float m = mean(data, length);
	
	for(int i = 0; i < length; i++) {
		diff = data[i] - m;
		sum += diff * diff;
	}
	return sum;
}

int remove_dups(float* data, int length)
{
	int size = length;
	for (int i = 0; i < size; i ++) {  
		for (int j = i + 1; j < size; j++){  
			if (data[i] == data[j])  
			{  
				for (int k = j; k < size - 1; k++)  
				{  
					data[k] = data[k + 1];  
				}  
				size--;  
				j--;      
			}  
		}  
	}  
	return size;	
}

int remove_zeros(float* data, int length)
{
	int size = length;
	int zeros = 0;
	for (int i = 0; i < size; i ++) {  
		if(data[i] != 0.0) {
			data[zeros++] = data[i]; 
		}
	}
	for(int i = zeros; i < size; i++) {
		data[i] = 0;
	}
	zeros = abs(zeros - size);
	return (size - zeros);	
}
