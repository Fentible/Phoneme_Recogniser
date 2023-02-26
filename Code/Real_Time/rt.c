#include "rt.h"

pthread_mutex_t lock;
pthread_cond_t cond;

short in_array[64000] = {0};
short pr_array[64000] = {0};

int in_size = 0;
int pr_size = 0;

void* process_thread(void* argp);
void* record_thread(void* argp);
void update_arrays(int in, int pr);
void boundary(void);

int STOP = 0;

short* read_wav(FILE* fp);

int main(void)
{
	dtw_init();
	pthread_mutex_init(&lock, NULL);
	pthread_cond_init(&cond, NULL);
		
	pthread_t tid;

	pthread_create(&tid, NULL, record_thread, NULL);
	pthread_create(&tid, NULL, process_thread, NULL);

	pthread_exit(NULL);
	while(1) {}
	
	return 0;
}

void* record_thread(void* argp)
{
	
	while(1) {
		sleep(1);
		record();
	}

	return NULL;
}

void* process_thread(void* argp)
{
	
	while(1) {
		// pr_count++;
		// check_in();
		sleep(1);
		update_arrays(1, 0);
	}

	return NULL;
}

/* void check_in(void) */
/* { */
/* 	if(in_size > 16) { */
/* 		// update_arrays(1, 0); */
/* 	} */
/* 	return; */
/* } */

void update_arrays(int in, int pr)
{
	pthread_mutex_lock(&lock); // lock the arrays
	if((pr == 0) && (in_size > 16)) {
		// move the in_buffer to the end of the process_buffer
		for(int i = 0; i < in_size; i++) {
			pr_array[i + pr_size] = in_array[i];
			
		}
		// zero out the entire buffer
		for(int i = 0; i < in_size; i++) {
			in_array[i] = 0;
		}
		// reset the size to allow appending to the start
		pr_size += in_size;
		in_size = 0;
		//for(int i = 0; i < pr_size; i++) {
		//	printf("pr_size[%d] = %d\n", i, pr_array[i]);
		// }
		if(pr_size > 64000 || in_size > 64000) {
			perror("Over 'pr' limit\n");
			exit(-1);
		}
		boundary(); //  { do_seperation -> dtw }
	} else {
		in_array[in_size] = in;
		in_size++;
		if(in_size > 64000 || pr_size > 64000) {
			perror("Over 'in' limit\n");
			exit(-1);
		}
	}
	pthread_cond_signal(&cond);
	pthread_mutex_unlock(&lock);

	return;
}

void record(void)
{
	DIR *p;
	struct dirent *pp;
	char* dir_name = "../Training/FEMALE/";

	if(!(p = opendir (dir_name))) {
		printf("Failed to open foder %s\n", dir_name);
		return;
	}
	char pathname[1024];
	FILE* fp;
	printf(":: Training folder :: %s\n", dir_name);
	if (p != NULL) {
		while ((pp = readdir (p)) != NULL) {
			
			short* sequence  = NULL;
			int length = strlen(pp->d_name);
			if (strncmp(pp->d_name + length - 4, ".wav", 4) == 0) {
				sprintf(pathname, "%s%s", dir_name, pp->d_name);
				fp = fopen(pathname, "rb");
				if(!fp) { printf("Error opening '%s' file", pathname); }
				printf("Opening '%s' file\n", pathname);
				sequence = read_wav(fp);
				int seq_length = sequence[0];
				for(int i = 1; i < seq_length; i++) {
					update_arrays(sequence[i], 1);
					usleep(62.5);
				}
				free(sequence);
				
			}
			
		}
		(void) closedir (p);
	}

	return;
}

short* read_wav(FILE* fp)
{

	char *buffer;
	long file_length;

	fseek(fp, 0, SEEK_END);
	file_length = ftell(fp);
	rewind(fp);
	buffer = (char *)malloc(file_length * sizeof(char));
	size_t r = fread(buffer, sizeof(char), file_length, fp);
	if(r != file_length) {
		perror("Fread fail");
		printf("Tried to read %ld but actually read %ld\n", file_length, r);
	}
	fclose(fp);

	int sample_bits =(int)((unsigned char)(buffer[35]) << 8 |
			       (unsigned char)(buffer[34]));

	int sample_size =(int)((unsigned char)(buffer[43]) << 24 |
			       (unsigned char)(buffer[42]) << 16 |
			       (unsigned char)(buffer[41]) << 8 |
			       (unsigned char)(buffer[40]));
	
	short* sequence = malloc((sample_size + 1) * sizeof(short));
	short byte = 0;
	if(sample_bits == 16) { // short size
		int incr = 44, pos = 0;
		while(incr != sample_size) {
			byte = (short)((unsigned char)(buffer[incr+1]) << 8|
					      (unsigned char)(buffer[incr]));
			memcpy(&sequence[pos], &byte, sizeof(byte));
			pos++;
			incr +=2;
		}
	}
	sequence[0] = sample_size;
	free(buffer);
	return sequence;

}

void boundary(void)
{
	int bound = 0, result = 0; // result = 0; // (send result to buzzer handler)
	while(1) {
		bound = next_boundary(pr_array, pr_size);
		if(bound == 0) {
			pr_size = shift_and_reduce(pr_array, pr_size, 32);
			pr_size = pr_size - 32;
			continue;
		}
		if(bound == -1) {
			break;
		}
		short* h = calloc(bound, sizeof(short));
		for(int i = 0; i < bound; i++) {
			h[i] = pr_array[i];
		}
		result = test_phoneme_utterance(h, bound);
		printf("Result :: %s\n", phones[result]->index->name);
		pr_size = shift_and_reduce(pr_array, pr_size, bound);
		pr_size = pr_size - bound;
	}
	

	return;
}
