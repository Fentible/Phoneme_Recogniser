#include "mfccs.h"

struct Phoneme** phones;

int truncation = 0;

void load_from_file(struct Phoneme* phone)
{
	char base[512];
	sprintf(base, "../device/%s/", phone->index->name);
	char length[512];
	long len = 0;
	char* end_ptr;
	DIR* theFolder = opendir(base);
	if(theFolder == NULL) {
		printf("Couldn't open folder :: %s\n", base);
		exit(-1);
	}
	struct dirent *next_file;
	while ( (next_file = readdir(theFolder)) != NULL )
	{
		if(strcmp(next_file->d_name, "..") == 0 || strcmp(next_file->d_name, ".") == 0
		   || strcmp(next_file->d_name, "device") == 0 || next_file->d_name[strlen(next_file->d_name) - 1] == 'f') {
			continue;
		}

		int n = 0, m = 0;
		// printf("File :: %s\n", next_file->d_name);
		while(next_file->d_name[n] != '_') {
			length[n] = next_file->d_name[n];
			n++;
		}
		length[n] = '\0';
		len = strtol(length, &end_ptr, 10);
		m = n + 1;
		while(next_file->d_name[m] != '\0') {
			length[m] = next_file->d_name[m];
			m++;
		}
		int exists = 0;
		for(int i = 0; i < phone->size_count; i++) {
			if(phone->size[i] == len) {
				phone->count[i]++;
				exists = 1;
			}
		}
		if(exists == 0) {
			phone->size = (int*)realloc(phone->size, (phone->size_count + 1) * sizeof(int));
			phone->count = (int*)realloc(phone->count, (phone->size_count + 1) * sizeof(int));
			phone->size[phone->size_count] = len;
			phone->count[phone->size_count] = 1;
			phone->size_count++;
			continue;
		}
	}
	closedir(theFolder);

	return;
}

void dtw_init(void)
{
	int j = 0;
	while(strcmp(p_codes[j], "\0") != 0) {
		j++;	
	}
	num_ph = j;

	phones = (struct Phoneme**)malloc(j * sizeof(struct Phoneme*));
	for(int i = 0; i < j; i++) {
		phones[i] = (struct Phoneme*)malloc(sizeof(struct Phoneme));
		phones[i]->size = (int*)malloc(sizeof(int));
		phones[i]->count = (int*)malloc(sizeof(int));
		phones[i]->mfcc = (float**)malloc(sizeof(float*));
		phones[i]->index = (struct Ph_index*)malloc(sizeof(struct Ph_index));
		phones[i]->index->i = i; 
		phones[i]->index->name = strdup(p_codes[i]);
		phones[i]->size_count = 0;
		phones[i]->use_count = 0;
		load_from_file(phones[i]);
	}
	 
	return;
}

struct Phoneme* get_mfcc(struct Phoneme* phone, int test_length)
{
	char base[128];
	// printf("Phone use count :: %d\n", phone->use_count);
	for(int i = 0; i < phone->use_count; i++) {
		free(phone->mfcc[i]);
	}
	phone->use_count = 0;
	
	sprintf(base, "../device/%s/", phone->index->name);
	// printf("Len :: %ld\n", len);
	char filename[512];

	for(int i = 0; i < phone->size_count; i++) {
		if(phone->size[i] != test_length)
			continue;
		for(int j = 0; j < phone->count[i]; j++) {
			sprintf(filename, "../device/%s/%d_%d.phn", phone->index->name, phone->size[i], j);
			phone->mfcc = (float**)realloc(phone->mfcc, (phone->use_count + 1) * sizeof(float*));
			phone->mfcc[phone->use_count] = (float*)calloc((phone->size[i] * floor(glbl_banks * glbl_test_trunc)), sizeof(float));
			FILE* fp = fopen(filename, "rb");
			if(fp == NULL) {
				printf("Failed to open :: %s\n", filename);
				continue;
			}
			char* end_ptr;
			char line[512];
			char *p; 
			int n = 0;
			while(fgets(line, sizeof(line), fp)) {
				p = strtok(line, " ");
				while((p = strtok(NULL, " ")) != NULL) {
					if(strcmp(p, "\n") == 0)
						continue;
					// printf("p :: '%s'\n", p);
					phone->mfcc[phone->use_count][n] = strtof(p, &end_ptr);
					n++;
					
				}
			}
			n = 0;
			phone->use_count++;
			fclose(fp);
		}
	}

	return phone;
}

/* struct Phoneme* get_mfcc(struct Phoneme* phone, int test_length) */
/* { */
/* 	char base[128]; */
/* 	// printf("Phone use count :: %d\n", phone->use_count); */
/* 	for(int i = 0; i < phone->use_count; i++) { */
/* 		free(phone->mfcc[i]); */
/* 	} */
/* 	phone->use_count = 0; */
	
/* 	sprintf(base, "../device/%s/", phone->index->name); */
/* 	char length[5]; */
/* 	long len = 0; */
/* 	char* end_ptr; */
/* 	DIR* theFolder = opendir(base); */
/* 	if(theFolder == NULL) { */
/* 		printf("Couldn't open folder :: %s\n", base); */
/* 		exit(-1); */
/* 	} */
/* 	struct dirent *next_file; */
/* 	while ( (next_file = readdir(theFolder)) != NULL ) */
/* 	{ */
/* 		if(strcmp(next_file->d_name, "../device/") == 0) */
/* 			continue; */

/* 		char filename[512]; */
/* 		sprintf(filename, "../device/%s/%s", phone->index->name, next_file->d_name); */
/* 		int n = 0, m = 0; */
/* 		while(next_file->d_name[n] != '_') { */
/* 			length[n] = next_file->d_name[n]; */
/* 			n++; */
/* 		} */
/* 		length[n] = '\0'; */
/* 		// printf("filename :: %s || length :: %s || test_length :: %d\n", next_file->d_name, length, test_length); */
/* 		len = strtol(length, &end_ptr, 10); */
/* 		m = n + 1; */
/* 		while(next_file->d_name[m] != '\0') { */
/* 			length[m] = next_file->d_name[m]; */
/* 			m++; */
/* 		} */
		
/* 		if(abs(len - test_length) <= 5) { */
/* 			// printf("Len :: %ld\n", len); */
/* 			phone->mfcc = (float**)realloc(phone->mfcc, (phone->use_count + 1) * sizeof(float*)); */
/* 			phone->mfcc[phone->use_count] = (float*)calloc(len, sizeof(float)); */
/* 			FILE* fp = fopen(filename, "rb"); */
/* 			char* end_ptr; */
/* 			char line[512]; */
/* 			char *p;  */
/* 			int n = 0; */
/* 			while(fgets(line, sizeof(line), fp)) { */
/* 				while( (p = strtok(line, " ")) != NULL) { */
/* 					phone->mfcc[phone->use_count][n] = strtof(p, &end_ptr); */
/* 					n++; */
/* 				} */
/* 			} */
/* 			phone->use_count++; */
/* 		} */
/* 	} */
/* 	closedir(theFolder); */
	
/* 	return phone; */
/* } */
