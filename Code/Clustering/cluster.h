#ifndef CLUSTER_H
#define CLUSTER_H

#include <memory.h>
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <unistd.h>

#include "../Misc/realloc.h"
#include "../Dynamic_Time_Warping/dtw.h"

struct Cluster* clust(float* data, int length, int amount, int offset, int num);
float SDCM(float* data, int length);

/** 
 * \struct Cluster
 * \brief Cluster structure :: used to hold the generated 1D Jenks Natural Breaks cluster
 * @centroids can be one or more centroids, these are generated from the @clust function
 * @values these are the values from the original data after it has been sorted and duplicates have been removed
 * the i'th array is associated with the i'th centroid in @centroids
 * @sizes is the length of each @values array with the corresponding index
 * @count is the number of centroids in @centroids
 * @gfv is the Goodness of Variance Fit and is a descriptor of how well the centroids describe the dataset
 */
struct Cluster {
	float* centroids;
	float** values;
	int* sizes;
	int count;
	float gfv;
} Cluster;

#endif
