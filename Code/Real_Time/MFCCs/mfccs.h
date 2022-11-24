#ifndef MFCCS_H
#define MFCCS_H

#include "../Misc/includes.h"

struct Phoneme* get_mfcc(struct Phoneme* phone, int length);
void dtw_init(void);

extern int truncation;
extern struct Phoneme** phones;
extern char* p_codes[];

#endif

