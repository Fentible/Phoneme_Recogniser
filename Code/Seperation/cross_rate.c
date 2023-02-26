/**
 * @file   cross_rate.c
 * @author T. Buckingham
 * @date   Tue May 17 00:23:35 2022
 * 
 * @brief  Zero crossing rate functions, primarily used in boundary detection and voice classification
 * 
 * 
 */
#include "cross_rate.h"

int is_positive(short num);
int f_is_positive(float num);

/** 
 * @brief Calculates the crossing rate of the given signal
 * 
 * @param signal The signal to process
 * @param signal_length The length of @signal
 * 
 * @return The total crossing rate of the signal
 */
int cross_rate(short* signal, int signal_length)
{

	int result = 0;
	for(int i = 0; i < signal_length - 1; i++) {
		if(is_positive(signal[i]) != is_positive(signal[i+1])) {
			result++;
			
		}
	}
	
	return result;
}

/** 
 * @brief Calculates the average cross rate of the given signal
 * 
 * @param signal The signal to process
 * @param signal_length The length of @signal
 * 
 * @return The averaged crossing rate as a floating point
 */
float favg_cross_rate(short* signal, int signal_length)
{

	float result = 0;
	for(int i = 0; i < signal_length - 1; i++) {
		if(is_positive(signal[i]) != is_positive(signal[i+1])) {
			result++;
		}
	}
	if(result == 0) 
		return result;
	else
		return (result / (signal_length - 1));
}

/** 
 * @brief Check if a short is postivie
 * 
 * @param num The number to check
 * 
 * @return 1 if positive, 0 if negative or zero.
 */
int is_positive(short num)
{
	if(num > 0) return 1;
	else return 0;
	
	return -1;
}

/** 
 * @brief The same as \fn cross_rate() but for float arrays as input
 * 
 * @param signal The signal to be processed
 * @param signal_length The length of @signal
 * 
 * @return The total cross rate of @signal
 */
float f_cross_rate(float* signal, int signal_length) 
{

	int result = 0;
	for(int i = 0; i < signal_length - 1; i++) {
		if(f_is_positive(signal[i]) != f_is_positive(signal[i+1])) {
			result++;
		}
	}

	return result;
}

/** 
 * @brief The same as \fn is_positive() but for floats
 * 
 * @param num The number to check if positive or not
 * 
 * @return 1 if positive, 0 if negative or 0
 */
int f_is_positive(float num)
{
	if(num > 0) return 1;
	else return 0;
	
	return -1;
}

/** 
 * @brief Calculates the short time zero cross average of a signal signal
 * 
 * @param signal The signal to process
 * @param signal_length The length of @signal
 * 
 * @return The short time averaged crossing rate as a floating point
 *
 * The signal is windowed using a Hanning window, the result from each window is summed then
 * the average of these windows is taken. Unlike \fn stavg_cross_rate() these windows do not overlap.
 */
float stavg_cross_rate_no_overlap(float* signal, int signal_length)
{
	float result = 0;

	int last = floor(signal_length / glbl_window_width);
        float** chunks =  hanning_chunks_no_overlap(signal, signal_length, glbl_window_width);

	double sum = 0, mean = 0, std = 0;
	
	for(int i = 0; i < last; i++) {
		result += f_cross_rate(chunks[i], glbl_window_width);	    
	}

	mean = result;
	for(int i = 0; i < last; i++) {
		sum += pow( ( (f_cross_rate(chunks[i], glbl_window_width)) - mean), 2);
	}
	if(sum == 0) {
		for(int i = 0; i < last; i++) {
			free(chunks[i]);
		}
		free(chunks);
		return 0;
	}
	std = sqrt(sum / (last + 1));
	result = 0;
	int tot = 0;
	for(int i = 0; i < last; i++) {
		double Z_score = ((f_cross_rate(chunks[i], glbl_window_width)) - mean) / std;
		if(Z_score < 2.0 && Z_score > -2.0) {
			result += f_cross_rate(chunks[i], glbl_window_width);
			tot++;
		}
	}
	result /= tot;
	
	for(int i = 0; i < last; i++) {
		free(chunks[i]);
	}
	free(chunks);
	return result;
}


/** 
 * @brief Calculates the short time zero cross average of a signal signal
 * 
 * @param signal The signal to process
 * @param signal_length The length of @signal
 * 
 * @return The short time averaged crossing rate as a floating point
 *
 * The signal is windowed using a Hanning window, the result from each window is summed then
 * the average of these windows is taken.
 */
float stavg_cross_rate(float* signal, int signal_length)
{
	float result = 0;
	int last = floor((signal_length - (glbl_window_width)) / (glbl_window_width / glbl_interval_div));
	
	float** chunks = hanning_chunks(signal, signal_length, glbl_window_width);

	for(int i = 0; i < last - 1; i++) {
		if(result == 0) {
			result = f_cross_rate(chunks[i], glbl_window_width);
		} else {
			result += f_cross_rate(chunks[i], glbl_window_width);
		}     
	}
	result = result / last;
	for(int i = 0; i < last; i++) {
		free(chunks[i]);
	}
	free(chunks);
	return result;
}
