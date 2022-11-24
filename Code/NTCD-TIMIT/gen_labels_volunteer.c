#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int main(void)
{

	FILE* fp = fopen("volunteer_labelfiles.mlf", "r");
	if(fp == NULL) {
		printf("Failed to open\n");
	}
	char line[256];
	char filename[512];
	char wav_filename[512];
	char outline[512];
	char subdir[512];
	char command[1024];
	char num[5];
	char* end_ptr;
	char *p;
	char *pp;
	long start = 0, end = 0;
	FILE* out;
	struct stat st = {0};
	if (stat("./labels/TRAIN", &st) == -1) {
		mkdir("./labels/TRAIN", 0700);
	}
	if (stat("./labels/TEST", &st) == -1) {
		mkdir("./labels/TEST", 0700);
	}
	while(fgets(line, sizeof(line), fp)) {
		//printf("line :: %s\n", line);
		p = strtok(line," ");
		if(p[0] == '.') {
			fclose(out);
			continue;
		}
		if(p[0] == '"') {
			// printf("%s\n", p);
			char gen = p[34];
			char num1 = p[32];
			char num2 = p[33];
			sprintf(num, "%c%c", num1, num2);
			long number = strtol(num, &end_ptr, 10);
			p = &p[54];
			p[strlen(p) - 6] = '\0';
			sprintf(wav_filename, "./clean/Clean/volunteers/%c%c%c/straightcam/%s.wav", num1, num2, gen, p);
			if(number < 45) {
				if(gen == 'M') {
					sprintf(subdir, "./labels/TRAIN/MALE/");
					struct stat st = {0};
					if (stat(subdir, &st) == -1) {
						mkdir(subdir, 0700);
					}
					sprintf(filename, "./labels/TRAIN/MALE/%s.PHN", p);
				} else {
					sprintf(subdir, "./labels/TRAIN/FEMALE/");
					struct stat st = {0};
					if (stat(subdir, &st) == -1) {
						mkdir(subdir, 0700);
					}
					sprintf(filename, "./labels/TRAIN/FEMALE/%s.PHN", p);
				}
			} else {
				if(gen == 'M') {
					sprintf(subdir, "./labels/TEST/MALE/");
					struct stat st = {0};
					if (stat(subdir, &st) == -1) {
						mkdir(subdir, 0700);
					}
					sprintf(filename, "./labels/TEST/MALE/%s.PHN", p);
				} else {
					sprintf(subdir, "./labels/TEST/FEMALE/");
					struct stat st = {0};
					if (stat(subdir, &st) == -1) {
						mkdir(subdir, 0700);
					}
					sprintf(filename, "./labels/TEST/FEMALE/%s.PHN", p);
				}
			}
			out = fopen(filename, "w");
			continue;
		}

		start = strtol(p, &end_ptr, 10);
	      	p = strtok(NULL, " ");
		end = strtol(p, &end_ptr, 10);
		p = strtok(NULL, " ");
		p[strlen(p) - 1] = '\0';
		// printf("%ld %ld %s\n", start, end, p);
		start = (start / 625);
		end = (end / 625);
		// printf("%ld %ld %s :: %s\n", start, end, p, filename);
		printf("Exporting :: %s\n", filename);
		printf("Moving %s to %s\n", wav_filename, subdir);
		sprintf(command, "cp %s %s", wav_filename, subdir);
		system(command);
		sprintf(outline, "%ld %ld %s\n", start, end, p);
		// printf("%s to %s\n", outline, filename);
		fprintf(out, "%s", outline);
	}



	return 0;
}
