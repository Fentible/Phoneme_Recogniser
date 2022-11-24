#!/bin/bash
# 95 : 750000 : 440 : -100 : 1250
#		   ( (change == -100 && ste_change > 1000 && entropy_change > 340) ||
#		     (change >= 200 && ste_change > 1500 && entropy_change > 400))) {

for(( i=4; i<=4; i=(i+1) )) # zc_incr
do
    for(( j=-80; j>=-90; j=(j-5) )) # ste_incr
    do
	for(( t=-85; t>=-95; t=(t-5) )) # entr_incr
	do
	    for(( n=425; n<=475; n=(n+10) )) # entr_incr
	    do
		for(( d=1050; d<=1450; d=(d+50) )) # ste_incr
		do

		    eval "find ./fails/ -name \"*.txt\" -type f -delete"
		    eval "find ./correct/ -name \"*.txt\" -type f -delete"
		    eval "find ./group/ -name \"*.txt\" -type f -delete"

		    eval "find ../Testing/TEST/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/MALE/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/FEMALE/res/ -name \"*.res\" -type f -delete"

		    eval "find ../Testing/SPKR1/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/SPKR1_NOSIL/res/ -name \"*.res\" -type f -delete"

		    eval "find ../Testing/SPKR1_CAFE/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/SPKR1_STRT/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/SPKR1_CAR/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/SPKR1_WHITE/res/ -name \"*.res\" -type f -delete"

		    eval "find ../MFCCs/ -name \"*.txt\" -type f -delete"
		    eval "./dtw.exe paa 2 window 128 banks 16 frame_limit 15 paa_op 0 dtw_window 200 interval_div 2 nfft 512 trunc 14 mfccs 99999 zc_incr ${i} ste_incr ${j} entr_incr ${t} larg_entr_incr ${n} larg_ste_incr ${d} test_iter 1 SIMPLE KNN knn 7 group_k 7 voice_k 7 ZC SPKR1 BOUNDS;"
		done
	    done
	done
    done
done






