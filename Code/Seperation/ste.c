/**
 * @file   ste.c
 * @author T. Buckingham
 * @date   Thu May  5 17:59:44 2022
 * 
 * @brief  Functions for calculating the short time energy of a signal
 * 
 * 
 */
#include "ste.h"

/** 
 * @brief Calculates the short time energy of a signal
 * 
 * @param signal The signal to process
 * @param signal_length The length of @param signal
 * @param window_length The length of the window to use
 * 
 * @return The calculated short time energy
 */
float short_time_energy(short* signal, int signal_length, int window_length)
{
	float result = 0;
	float* window = hanning_window(window_length);
	float* ste_buff = calloc(signal_length, sizeof(float));
	int pos = 0, incr = 1;
	for(int i = 0; i < signal_length; i++) {
		ste_buff[i] = (signal[i] * signal[i]); // * window[pos];
		pos++;
		if(pos == glbl_window_width) {
			pos = 0;
			incr++;
		}
	}

	for(int i = 0; i < signal_length; i++) {
		result += ste_buff[i];
		
	}
	result /= incr;
	free(window);
	free(ste_buff);
	return result;
}

/** 
 * @brief Calculates the short time energy of a signal; this version is for
 * a float input signal
 * 
 * @param signal The signal to process
 * @param signal_length The length of @param signal
 * @param window_length The length of the window to use
 * 
 * @return The calculated short time energy
 */
float f_short_time_energy(float* signal, int signal_length, int window_length)
{
	float result = 0;
	float* window = hanning_window(window_length);
	float* ste_buff = calloc(signal_length, sizeof(float));
	int pos = 0, incr = 1;
	for(int i = 0; i < signal_length; i++) {
		ste_buff[i] = (signal[i] * signal[i]); // * window[pos];
		pos++;
		if(pos == glbl_window_width) {
			pos = 0;
			incr++;
		}
	}

	for(int i = 0; i < signal_length; i++) {
		result += ste_buff[i];
		
	}
	result /= incr;
	free(window);
	free(ste_buff);
	return result;
}


		
