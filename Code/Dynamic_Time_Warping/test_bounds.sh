#!/bin/bash
for(( i=93; i<=97; i=(i+10) )) # zc_incr
do
    for(( j=650000; j<=850000; j=(j+50000) )) # ste_incr
    do
	for(( t=380; t<=420; t=(t+10) )) # entr_incr
	do
	    for(( n=-90; n>=-100; n=(n-2) )) # neg_entr_incr
	    do
		for(( l=1050; l<=1250; l=(l+150) )) # larg_entr_incr
		do
		    
		    eval "find ./aao/fails/ -name \"*.txt\" -type f -delete"
		    eval "find ./aao/correct/ -name \"*.txt\" -type f -delete"
		    eval "find ./aao/group/ -name \"*.txt\" -type f -delete"

		    eval "find ./fails/ -name \"*.txt\" -type f -delete"
		    eval "find ./correct/ -name \"*.txt\" -type f -delete"
		    eval "find ./group/ -name \"*.txt\" -type f -delete"

		    eval "find ../Testing/TEST/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/TEST/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/TEST/res/ -name \"*.res\" -type f -delete"

		    eval "find ../Testing/MALE/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/MALE/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/MALE/res/ -name \"*.res\" -type f -delete"

		    eval "find ../Testing/FEMALE/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/FEMALE/res/ -name \"*.res\" -type f -delete"
		    eval "find ../Testing/FEMALE/res/ -name \"*.res\" -type f -delete"

		    eval "find ../MFCCs/ -name \"*.txt\" -type f -delete"
		    
		    eval "./dtw.exe paa 2 window 128 banks 40 paa_op 0 dtw_window 200 interval_div 2 nfft 512 trunc 24 mfccs 99999 SIMPLE zc_incr 93 ste_incr 650000 entr_incr 420 neg_incr -94 larg_incr 1200 FEMALE KNN knn 7 test_iter 1 GROUP ZC;"

		done
	    done
	done
    done
done





