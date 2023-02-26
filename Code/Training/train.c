/**
 * \file train.c
 * \Handles training of the phonemes and reading of .wav files
 */
#include "train.h"

struct Phoneme** phones; /* A program scope list of all phonemes using for comparions */

/**
 * @p_codes contains a list of phoneme codes which can be found within the given dataset
 * TODO allow the loading of this list from an external config file
 * This program is designed for the TCD and NTCD TIMIT data-sets and as such some values have
 * been hard coded. To test grouped and/or voiced KNN with a different set of phonemes the
 * values in @file knn.c need to be changed per function.
 */
char* p_codes[] = {"\0", "b", "d", "k", "p", "t", "g",                                                       // Stops [1 - 6] - OBSTRUENT
	                 "jh", "ch",                                                                         // Affri [7 - 8] - OBSTRUENT
			 "s", "sh", "th", "v", "f", "dh", "z",                                               // Frics [9 - 15] - OBSTRUENT
			 "m", "n", "ng",                                                                     // Nasal [16 - 18] - SONORANT
			 "l", "r", "hh", "w", "y",                                                           // Semiv [19 - 23] - SONORANT
	                 "aa", "ae", "ah", "aw", "er", "ay", "eh", "ey", "ih", "iy", "ow", "oy", "uh", "uw", // Vowel [24 - 37] - SONORANT
			 "sil",                                                                              // Other [38] - OTHER
	           "\0"};

/**
 * @p_group contains a list of phoneme groups 
 */
char* p_group[] = {"STOP", "AFRI", "FRIC", "NASL", "SEMV", "VOWL", "OTHR"};

/**
 * @OBSTR @SONOR @OTHER are types of voiceness a phoneme (or sil) can be
 */
int OBSTR = 0;
int SONOR = 1;
int OTHER = 2;

/**
 * @ph_zc_max @ph_zc_min @ste_min @ste_max
 * All were used to determine if a phoneme should be tested using DTW or KNN
 * if the test's values were not within the range of the phoneme's values then
 * no test would be done
 */ 
float* ph_zc_max; 
float* ph_zc_min;
float* ste_min;
float* ste_max;


int limit_changed = 0; /* Keeps track of how many limits were changed \deprecated as the window is now based on phoneme's length  */


int num_tt = 0; /* TODO what is this? */
int num_ph = 0; /* Stores the number of phonemes from @p_codes */

/**
 * @max_sil @min_sil are values that store the maximum and minimum raw signal values found in testing
 * \brief used when testing numerous silence detection methods
 */
int max_sil = 0;
int min_sil = INT_MAX;
/**
 * @max_sil_dB @min_sil_dB are values that store the maximum and minimum raw signal values in decibels found in testing
 * \brief used when testing numerous silence detection methods
 */
float max_sil_dB = 0;
float min_sil_dB = INT_MAX;
/**
 * @max_sil_zc @min_sil_zc are values that store the maximum and minimum zero cross values found in testing
 * \brief used when testing numerous silence detection methods
 */
float max_sil_zc = 0;
float min_sil_zc = INT_MAX;
/**
 * @max_sil_flt @min_sil_flt are values that store the maximum and minimum flatness values found in testing
 * \brief used when testing numerous silence detection methods
 */
float max_sil_flt = 0;
float min_sil_flt = INT_MAX;
/**
 * @max_sil_ste @min_sil_ste are values that store the maximum and minimum short time energy values found in testing
 * \brief used when testing numerous silence detection methods
 */
float max_sil_ste = 0;
float min_sil_ste = INT_MAX;
/**
 * @max_sil_mean @min_sil_mean are values that store the maximum and minimum mean values found in testing
 * \brief used when testing numerous silence detection methods
 */
float max_sil_mean = 0;
float min_sil_mean = INT_MAX; 

static char* male_dir = "../Training/MALE/";
static char* female_dir = "../Training/FEMALE/";
static char* train_dir = "../Training/TRAIN/";
static char* spkr1_dir = "../Training/SPKR1/";
static char* spkr1_nosil_dir = "../Training/SPKR1_NOSIL/";

void update_zc(float* signal, int signal_length, int n);
void update_ste(short* signal, int signal_length, int n);
void update_sil(short* signal, int signal_length);
void update_sil_db(short* signal, int signal_length);
void update_sil_mean(short* signal, int signal_length);
void update_sil_flat(short* signal, int signal_length, char* filename);
void update_sil_ste_frame(short* array, int length, char* filename);

/**
 * \fn mfcc_size()
 * \brief Returns the total number of values in an MFCC given the raw time domain signal
 * This function returns the length of the MFCC if the MFCC was a 1D array so if the number of
 * frames is needed then dividing the result by the number of coefficients or using \fn frame_amount() is needed
 */
int mfcc_size(int signal_length)
{

	int width = floor(signal_length / glbl_paa);
	int amount = floor((width - (glbl_window_width)) / (glbl_window_width / glbl_interval_div));
	if(amount <= 0) {
		amount = floor(width / glbl_window_width);
	}
	if(amount > glbl_frame_limit) {
		amount = glbl_frame_limit;
	}
	int new_size = amount * floor((glbl_banks) * glbl_test_trunc);
	return new_size;
	
}

/**
 * \fn frame_amount()
 * \brief Returns the number of frames an MFCC will have from a time domain signal of the given length
 */
int frame_amount(int signal_length)
{

	int width = floor(signal_length / glbl_paa);
	int amount = floor((width - (glbl_window_width)) / (glbl_window_width / glbl_interval_div));
	if(amount <= 0) {
		amount = floor(width / glbl_window_width);
	}
	if(amount > glbl_frame_limit) {
		amount = glbl_frame_limit;
	}
	return amount;	
	
}

/**
 * \fn cubic_interpolate()
 * \brief Interpolates between points @y1 and @y2 using @y0 and @y3 to apply the cubic function
 * @mu is the distance between @y1 and @y2 the user wishes the return to be and so a value between 0 and 1 is used
 */ 
double cubic_interpolate(short y0, short y1, short y2, short y3, double mu)
{
	short a0 = 0, a1 = 0, a2 = 0, a3 = 0, mu2  = 0;

	mu2 = mu*mu;
	a0 = y3 - y2 - y0 + y1;
	a1 = y0 - y1 - a0;
	a2 = y2 - y0;
	a3 = y1;
	return(a0*mu*mu2+a1*mu2+a2*mu+a3);
}

long double ld_cubic_interpolate(long double y0, long double y1, long double y2, long double y3, long double mu)
{
	long double a0 = 0, a1 = 0, a2 = 0, a3 = 0, mu2  = 0;

	mu2 = mu*mu;
	a0 = y3 - y2 - y0 + y1;
	a1 = y0 - y1 - a0;
	a2 = y2 - y0;
	a3 = y1;
	return(a0*mu*mu2+a1*mu2+a2*mu+a3);
}

/**
 * \fn train()
 * \brief is the main control function for training
 * \fn train() handles reading in all files from the specified folder, a .wav is found a read then a corresponding
 * .PHN file found. These are passed to \fn allocate_ph() to read, MFCC'd and applied to a phoneme prototype
 */
void train(void)
{
	
	DIR *p;
	struct dirent *pp;
	char dir_name[35];
	if(MALE) {
		memcpy(dir_name, male_dir, sizeof(char) * (strlen(male_dir) + 1));
	} else if (FEMALE) {
		memcpy(dir_name, female_dir, sizeof(char) * (strlen(female_dir) + 1));
	} else if (SPKR1) {
		memcpy(dir_name, spkr1_dir, sizeof(char) * (strlen(spkr1_dir) + 1));
	} else if (SPKR1_NOSIL) {
		memcpy(dir_name, spkr1_nosil_dir, sizeof(char) * (strlen(spkr1_nosil_dir) + 1));
	} else {
		memcpy(dir_name, train_dir, sizeof(char) * (strlen(train_dir) + 1));
	}
	if(!(p = opendir (dir_name))) {
		printf("Failed to open foder %s\n", dir_name);
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
	
	printf(":: Training folder :: %s :: %d files to train\n", dir_name, files_to_test);
	if (p != NULL) {
		while ((pp = readdir (p)) != NULL) {
			
			short* sequence  = NULL;
			struct wav_file* wv = NULL;
			int length = strlen(pp->d_name);
			if (strncmp(pp->d_name + length - 4, ".wav", 4) == 0 || strncmp(pp->d_name + length - 4, ".WAV", 4) == 0) {
				int len = 8;
				if(strncmp(pp->d_name + length - 4, ".wav", 4) == 0)
					len = 4;
				sprintf(pathname, "%s%s", dir_name, pp->d_name);
				fp = fopen(pathname, "rb");
				if(!fp) { printf("Error opening '%s' file", pathname); }
				wv = read_wav(fp);
				sequence = wv->seq;
			        pp->d_name[length-len] = '\0';
				strcat(pp->d_name, ".PHN");
				sprintf(pathname, "%s%s", dir_name, pp->d_name);
				fp = fopen(pathname, "r");
				if(!fp) { printf("Error opening '%s' file", pathname); }
				allocate_ph(fp, sequence, 0);
				free(sequence);
				
			}
			free(wv);			
		}
		(void) closedir (p);
	}
}

/**
 * \fn resize()
 * \brief resizes a given signal (@shorter) from size @s to size @l using cubic interpolation
 */ 
short* resize(short* shorter, size_t s, size_t l)
{
	short* new_sequence = calloc(l, sizeof(short));
	for(unsigned int i = 0; i < s; i++) {
		new_sequence[i] = shorter[i];
	}
	free(shorter);
	size_t i = floor(s / 4);
	while(s < l) {
		for(size_t j = s-1; j >= i ; j--) {
			new_sequence[ j + 1 ] = new_sequence[ j ];
		}
	        		
	        new_sequence[i] = cubic_interpolate((short)new_sequence[i - 2], (short)new_sequence[i - 1], (short)new_sequence[i + 1], (short)new_sequence[i + 2], 0.5);
		s++;
		if(i + 2 >= s - 5){  i = 5; }
		i+=2;	
	}
	return new_sequence;

	
}

long double* ld_resize(long double* shorter, size_t s, size_t l)
{
	long double* new_sequence = calloc(l, sizeof(long double));
	for(unsigned int i = 0; i < s; i++) {
		new_sequence[i] = shorter[i];
	}
	free(shorter);
	size_t i = floor(s / 2);
	while(s < l) {
		for(size_t j = s-1; j >= i ; j--) {
			new_sequence[ j + 1 ] = new_sequence[ j ];
		}
		
	        new_sequence[i] = ld_cubic_interpolate((long double)new_sequence[i - 2], (long double)new_sequence[i - 1], (long double)new_sequence[i + 1], (long double)new_sequence[i + 2], 0.5);
		s++;
		if(i + 2 >= s - 5){  i = 5; }
		i+=2;	
	}
	return new_sequence;
}

/**
 * \fn train_ph_mfcc()
 * \brief the given sequence is added to the given phoneme prototype
 * This function handles adding new signals to the phoneme prototype
 * the raw time domain signal is added to an array, the size is added and the count of
 * number of signals in incremented. If the number of signals is above the given limit @glbl_mfccs
 * then the signals are interpolated to the nearest and averaged
 * @new the new sequences length - will be compared using the \fn mfcc_length()
 * @phone which phoneme the signal should be added to
 * @sequence the raw time domain signal
 */
short* train_ph_mfcc(int new, short* sequence, struct Phoneme* phone)
{

	if(new == 0 || mfcc_size(new) == 0) {
		printf("Zero size not permitted :: new : %d || mfcc(new) : %d\n", new, mfcc_size(new));
		return sequence;
	}

	if(phone->size_count == 0) {
		phone->sequence = (long double**)calloc(phone->size_count + 1, sizeof(long double*));	
		phone->sequence[0] = (long double*)calloc(new, sizeof(long double));
		for(int i = 0; i < new; i++) {
			phone->sequence[0][i] = sequence[i];
		}	  
	// 	memcpy(phone->sequence[0], sequence, sizeof(short) * new);
		phone->size[0] = new;
		phone->amounts[0] = 1;
		phone->size_count++;
		return sequence;
	}
	for(int i = 0; i < phone->size_count; i++) { 
		if(mfcc_size(new) == mfcc_size(phone->size[i]) && phone->size_count >= glbl_mfcc_num) {
			if(new > phone->size[i]) {
				for(int j = 0; j < phone->size[i]; j++) {
					long double temp = phone->sequence[i][j];
					phone->sequence[i][j] += sequence[j];
					if(temp > phone->sequence[i][j] && sequence[j] > 0) {
						printf("Overflow in training! : %Lf -> %Lf\n", temp, phone->sequence[i][j]);
						exit(-1);
					}
					if(temp < phone->sequence[i][j] && sequence[j] < 0) {
						printf("Underflow in training! : %Lf -> %Lf\n", temp, phone->sequence[i][j]);
						exit(-1);
					}
				}
			} else if(phone->size[i] > new) {
				for(int j = 0; j < new; j++) {
					long double temp = phone->sequence[i][j];
					phone->sequence[i][j] += sequence[j];
					if(temp > phone->sequence[i][j] && sequence[j] > 0) {
						printf("Overflow in training! : %Lf -> %Lf\n", temp, phone->sequence[i][j]);
						exit(-1);
					}
					if(temp < phone->sequence[i][j] && sequence[j] < 0) {
						printf("Underflow in training! : %Lf -> %Lf\n", temp, phone->sequence[i][j]);
						exit(-1);
					}
				}
			} else {
				for(int j = 0; j < phone->size[i]; j++) {
					long double temp = phone->sequence[i][j];
					phone->sequence[i][j] += sequence[j];
					if(temp > phone->sequence[i][j] && sequence[j] > 0) {
						printf("Overflow in training! : %Lf -> %Lf\n", temp, phone->sequence[i][j]);
						exit(-1);
					}
					if(temp < phone->sequence[i][j] && sequence[j] < 0) {
						printf("Underflow in training! : %Lf -> %Lf\n", temp, phone->sequence[i][j]);
						exit(-1);
					}
				}
			}
		        phone->amounts[i]++;
			return sequence;
		} 
	}

	// increase size of sequence**
	// add new sequence*
	// add new size*
	// add and increment new amounts
	// increment size_count
	if(phone->size_count < glbl_mfcc_num) {
		phone->sequence = realloc(phone->sequence, (phone->size_count + 1) * sizeof(long double*));
		
		phone->sequence[phone->size_count] = (long double*)calloc(new, sizeof(long double));
		for(int i = 0; i < new; i++) {
		 	phone->sequence[phone->size_count][i] = sequence[i];
		}
		//	memcpy(phone->sequence[phone->size_count], sequence, sizeof(short) * new);
		
		phone->size = (int*)realloc(phone->size, sizeof(int) * (phone->size_count + 1));
		phone->size[phone->size_count] = new;
		
		phone->amounts = (int*)realloc(phone->amounts, sizeof(int) * (phone->size_count + 1));
		phone->amounts[phone->size_count] = 1;
		phone->size_count++;	

	} else {
		int smallest_diff = 0;
		int p = 0;
		for(int i = 0; i < phone->size_count; i++) {
			if(i == 0 || abs(new - phone->size[i]) < smallest_diff) {
				p = i;
				smallest_diff = abs(new - phone->size[i]);
			}

		}
		if(new > phone->size[p]) {
			phone->sequence[p] = ld_resize(phone->sequence[p], phone->size[p], new);
			for(int j = 0; j < phone->size[p]; j++) {
				phone->sequence[p][j] += sequence[j];
			}
		} else if(phone->size[p] > new) {
			sequence = resize(sequence, new, phone->size[p]);
			for(int j = 0; j < new; j++) {
				phone->sequence[p][j] += sequence[j];
			}
		} else {
			for(int j = 0; j < phone->size[p]; j++) {
				phone->sequence[p][j] += sequence[j];
			}

		}
		phone->amounts[p]++;
		
	}
		
	return sequence;
}

/**
 * \fn allocate_ph()
 * \brief This function is used by \fn train() and \fn test() to read phoneme sequences from .wav files using .PHN data
 * If training then the signal is passed to \fn train_ph_mfcc() and min max values for silence, ste and zc values are updated here as well
 * If testing the signal for the phoneme is found then passed to \fn test_phoneme()
 */
void allocate_ph(FILE* fp, short* wav, unsigned char t_t)
{
	char line[256];
	char* end_ptr;
	char *p; 
	long start_end[2];
	
	while(fgets(line, sizeof(line), fp)) {
		p = strtok(line," ");
		start_end[0] = strtol(p, &end_ptr, 10);
	      	p = strtok(NULL, " ");
		start_end[1] = strtol(p, &end_ptr, 10);
		p = strtok(NULL, " ");
		while(p[strlen(p) - 1] == '\n' || p[strlen(p) - 1] == ' ') {
			p[strlen(p) - 1] = '\0';
		}

		if((start_end[1] - start_end[0]) <= 0) {
			printf("Negative length from wav %s :: %ld - %ld\n", p, start_end[0], start_end[1]);
			continue;
		}
		if(strcmp(p, "sil") == 0) {
			start_end[1] = (start_end[1] + start_end[0]) / 2;
		}
		short* h = (short*)calloc((start_end[1] - start_end[0]), sizeof(short));
		for(int n = 1; n < num_ph ; n++) {
			if(strcmp(p_codes[n], p) == 0) {
				int m  = 0;
				// printf("length :: %ld\nstart :: %ld\nend :: %ld\nph :: %s\n", (start_end[1] - start_end[0]), start_end[0], start_end[1], p);
				for(int j = start_end[0]; j < start_end[1]; j++) {
					h[m] = wav[j];
					m++;	
				}
				int length = start_end[1] - start_end[0];
				
				if(length <= 0) {
					printf("Negative length from wav :: %ld - %ld\n", start_end[0], start_end[1]);
					break;
				}
				for(int i = 1; i < length; i++) {
					h[i] = h[i] - (0.95 * h[i-1]);
				}
				int signal_length = start_end[1] - start_end[0] - 1;
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
				// h = resize(h, signal_length, (signal_length * 2));
				// signal_length = (signal_length * 2);
				if(t_t == TRAIN) {
					// int end =  floor((start_end[1] - start_end[0] - 1));
					
					// update_pca(h, signal_length, n); 
					float* signal = (float*)calloc(signal_length, sizeof(float));
					for(int i = 0; i < signal_length; i++) {
						signal[i] = h[i];
					}
					int MAXSIZE = 0xFFF;
					char proclnk[0xFFF];
					char filename[0xFFF];
					int fno;
					ssize_t r;
					
					fno = fileno(fp);
					sprintf(proclnk, "/proc/self/fd/%d", fno);
					r = readlink(proclnk, filename, MAXSIZE);
					filename[r] = '\0';
					
					update_zc(signal, signal_length, n);
					update_ste(h, signal_length, n);
					
					if(strcmp("sil", phones[n]->index->name) == 0) {
						update_sil_zc(h, signal_length, filename);
						update_sil_ste(h, signal_length, filename);
					}
					// printf("fp -> fno -> filename: %p -> %d -> %s\n",
					// (void*)fp, fno, filename);
					h = train_ph_mfcc(signal_length, h, phones[n]);
					trained++;
					phones[n]->trained++;
										
					if(trained % 5000 == 0)
						printf("::       TRAINED P#       ::  %02d:%02d:%02d  ::  %05d\n", hour(difftime(time(NULL), start)), minu(difftime(time(NULL), start)), seco(difftime(time(NULL), start)), trained);
					free(signal);
				} else {
					per_correct[n]++;
					test_phoneme(h, signal_length, p);
					// test_phoneme_pca(h, start_end, p);
					// test_phoneme_aao(h, start_end, p);
				}
			}
		}
		if(h != NULL) {
			free(h);
		}
		
	}

	fclose(fp);

	return;
}

/**
 * \fn read_wav()
 * \brief reads the wav file in full and returns a \struct wav_file containing the data and length of data
 */
struct wav_file* read_wav(FILE* fp)
{

	char *buffer;
	long file_length;
	struct wav_file* wv = (struct wav_file*)malloc(sizeof(struct wav_file));
	fseek(fp, 0, SEEK_END);
	file_length = ftell(fp);
	rewind(fp);
	buffer = (char *)malloc(file_length * sizeof(char));
	fread(buffer, file_length, 1, fp);
	fclose(fp);

	int sample_bits =(int)((unsigned char)(buffer[35]) << 8 |
			       (unsigned char)(buffer[34]));

	int sample_size =(int)((unsigned char)(buffer[43]) << 24 |
			       (unsigned char)(buffer[42]) << 16 |
			       (unsigned char)(buffer[41]) << 8 |
			       (unsigned char)(buffer[40]));
	
	short* sequence = calloc((sample_size / 2), sizeof(short));
	short byte = 0;
	if(sample_bits == 16) { // short size
		int incr = 44, pos = 0;
		while(incr != (sample_size + 44)) {
			byte = (short)((unsigned char)(buffer[incr+1]) << 8|
					      (unsigned char)(buffer[incr]));
			memcpy(&sequence[pos], &byte, sizeof(byte));
			pos++;
			incr +=2;
		}
	} else {
		printf("Samples bits != 16...\n");
		exit(-1);
	}
	wv->seq = sequence;
	wv->length = (sample_size / 2);
	free(buffer);
	return wv;

}

void update_sil_mean(short* signal, int signal_length)
{
	double total = 0, mean = 0;
	for(int i = 0; i < signal_length - 2; i++) {
		total += signal[i];
	}
	mean = total / (signal_length - 2);
	
	for(int i = 0; i < signal_length - 2; i++) {
		// double Z_score = (signal[i] - mean) / std;
		if(mean < min_sil_mean) {
			min_sil_mean = mean;
		}
		if(mean > max_sil_mean) {
			max_sil_mean = mean;
		}
	}

	return;
}

void update_sil_ste(short* array, int length, char* filename)
{
	float* signal = (float*)calloc(length, sizeof(float));
	for(int i = 0; i < length; i++) {
		signal[i] = array[i];
	}
	float ste = f_short_time_energy(signal, length, glbl_window_width);
	FILE* fp = fopen("./ste_.txt", "a");
	if(ste > max_sil_ste) 
		max_sil_ste = ste;
	if(ste < min_sil_ste)
		min_sil_ste = ste;

	// fprintf(fp, "%f : %s\n", ste, filename);
	free(signal);
	fclose(fp);
	return;
}

void update_sil_zc(short* array, int length, char* filename)
{
	float n = 0;
	float zc = 0;
	FILE* fp = fopen("./zc_.txt", "a");
	for(int i = 0; i < length - glbl_window_width; i+=glbl_window_width) {
		zc += favg_cross_rate(&array[i], glbl_window_width);
		n++;
	}
	zc /= n;
	if(zc > max_sil_zc) 
		max_sil_zc = 0.30;
	if(zc < min_sil_zc)
		min_sil_zc = 0.22;
	
	// fprintf(fp, "%f : %s\n", zc, filename);
	fclose(fp);
	return;
}

void update_sil(short* signal, int signal_length)
{
	double sum = 0, total = 0, mean = 0, std = 0;
	for(int i = 0; i < signal_length - 2; i++) {
		total += signal[i];
	}
	mean = total / (signal_length - 2);
	for(int i = 0; i < signal_length - 2; i++) {
		sum += pow( (signal[i] - mean), 2);
	}
	std = sqrt(sum / (signal_length - 2));
	for(int i = 0; i < signal_length - 2; i++) {
		double Z_score = (signal[i] - mean) / std;		
		if(signal[i] < min_sil && Z_score < 0.1 && Z_score > -0.1) {
			min_sil = signal[i];
		}
		if(signal[i] > max_sil && Z_score < 0.1 && Z_score > -0.1) {
			max_sil = signal[i];
		}
	}

	return;
}

float flatness(float* chunk, int length)
{
	float result = 0, sum = 0, log_sum = 0;
	float N = 1.0/length;
	for(int i = 0; i < length; i++) {
		sum += chunk[i];
		if(log_sum != 0)
			log_sum += log(chunk[i]);
	}
	if(N == 0 || sum == 0) {
		return 0;
	}
	result = exp(N * log_sum) / (N * sum);
	return result;
}

void update_sil_flat(short* signal, int signal_length, char* filename)
{
	FILE* fp = fopen("./flat_.txt", "a");
	float* sig = (float*)calloc(signal_length, sizeof(float));
	for(int i = 0; i < signal_length; i++) {
		sig[i] = signal[i];
	}
	float flat = flatness(sig, signal_length);
	if(flat < min_sil_flt)
		min_sil_flt = 0;
	if(flat > max_sil_flt)
		max_sil_flt = 1.1;

	// fprintf(fp, "%f : %s\n", flat, filename);
	free(sig);
	fclose(fp);
	return;
}

void update_sil_db(short* signal, int signal_length)
{
	float sum = 0, total = 0, mean = 0, std = 0;
	for(int i = 0; i < signal_length - 2; i++) {
		total += signal[i];
	}
	mean = total / (signal_length - 2);
	for(int i = 0; i < signal_length - 2; i++) {
		sum += pow( (signal[i] - mean), 2);
	}
	std = sqrt(sum / (signal_length - 2));
	for(int i = 0; i < signal_length - 2; i++) {
		float Z_score = (signal[i] - mean) / std;
		float amplitude = abs(signal[i]) / 32767;
		if(amplitude == 0)
			amplitude = FLT_EPSILON;
		float dB = 20 * log10(amplitude);
		if(dB < min_sil_dB && Z_score < 3 && Z_score > -3) {
			min_sil_dB = 10;
		}
		if(dB > max_sil_dB && Z_score < 3 && Z_score > -3) {
			max_sil_dB = -10;
		}
	}

	return;
}

void update_zc(float* signal, int signal_length, int n)
{
	
	float stavgnoc = stavg_cross_rate_no_overlap(signal, signal_length);
	if(fpclassify(stavgnoc) == FP_INFINITE || fpclassify(stavgnoc) == FP_NAN) {
		printf("Stavgnoc is not usable :: %.3f\n", stavgnoc);
	}
		
	if(phones[n]->trained == 0 || stavgnoc > ph_zc_max[n]) {
		ph_zc_max[n] = stavgnoc;
	}
	if(phones[n]->trained == 0 || stavgnoc < ph_zc_min[n]) {
		ph_zc_min[n] = stavgnoc;
	}

	return;
}

void update_ste(short* signal, int signal_length, int n)
{

	double ste = short_time_energy(signal, signal_length, glbl_window_width);
	if(fpclassify(ste) == FP_INFINITE || fpclassify(ste) == FP_NAN) {
		printf("Ste is not usable :: %.3f\n", ste);
	}

	if(phones[n]->trained == 0 || ste > ste_max[n]) {
		ste_max[n] = ste;
	}
	if(phones[n]->trained == 0 || ste < ste_min[n]) {
		ste_min[n] = ste;
	}	

	return;
}

void update_sil_ste_frame(short* array, int length, char* filename)
{
	FILE* fp = fopen("./ste_frame.txt", "a");
	float ste = 0;
	float* signal = (float*)calloc(length, sizeof(float));
	for(int i = 0; i < length; i++) {
		signal[i] = array[i];
	}
	for(int i = 0; i < (length - glbl_window_width); i+=glbl_window_width) {
		ste = f_short_time_energy(&signal[i], glbl_window_width, glbl_window_width);
		if(ste > max_sil_ste && (max_sil_ste * 1.25) > ste)
			max_sil_ste = ste;
			
	}
	fprintf(fp, "%f : %s\n", ste, filename);
	fclose(fp);
	return;
}

	      
