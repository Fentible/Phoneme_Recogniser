#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

int incr = 0;

int main(void)
{

	FILE* fp = fopen("lipspeaker_labelfiles.mlf", "r");
	if(fp == NULL) {
		printf("Failed to open\n");
	}
	char line[256];
	char filename[512];
	char outline[512];
	char subdir[512];
	char num[5];
	char* end_ptr;
	char *p; 
	long start = 0, end = 0;
	FILE* out;
	while(fgets(line, sizeof(line), fp)) {
		printf("line :: %s\n", line);
		p = strtok(line," ");
		if(p[0] == '.') {
			continue;
		}
		if(p[0] == '"') {
			if(p[40] != '1')
				continue;
			printf("%s\n", p);
			char gen = p[34];
			char num1 = p[32];
			char num2 = p[33];
			sprintf(num, "%c%c", num1, num2);
			int fold = atoi(num);
			p = &p[54];
			p[strlen(p) - 6] = '\0';
			incr++;
			if(incr >= 113) {
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
		sprintf(outline, "%ld %ld %s\n", start, end, p);
		printf("%s to %s\n", outline, filename);
		fprintf(out, "%s", outline);
	}



	return 0;
}
