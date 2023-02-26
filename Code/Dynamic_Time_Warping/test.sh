#!/bin/bash
for i in 128 # Hanning window width
do 
    for j in 16 # Number of filter banks
    do
	m=$(($j/2))
	mnm=$(($j-1))
	for((t=14;t<=14;t=t+6)) # Number of coefficients
	do
	    for ((n=200;n<=200;n+=200)) # DTW limit
	    do
		for p in 512 # NFFT
		do
		    for d in 2 # Window overlap division
		    do
			for k in 15 # Base KNN values
			do
			    for v in 15 # Frame limit
			    do
				for a in 2
				do


				    eval "find ./fails/ -name \"*.txt\" -type f -delete"
				    eval "find ./correct/ -name \"*.txt\" -type f -delete"
				    eval "find ./group/ -name \"*.txt\" -type f -delete"
				    
				    eval "find ../Testing/TEST/res/ -name \"*.res\" -type f -delete"
				    eval "find ../Testing/MALE/res/ -name \"*.res\" -type f -delete"
				    eval "find ../Testing/FEMALE/res/ -name \"*.res\" -type f -delete"
				    # eval "find ../Testing/SPKR1/res/ -name \"*.res\" -type f -delete"
				    # eval "find ../Testing/SPKR1_NOSIL/res/ -name \"*.res\" -type f -delete"
				    # eval "find ../Testing/SPKR1_CAR/res/ -name \"*.res\" -type f -delete"
				    # eval "find ../Testing/SPKR1_STRT/res/ -name \"*.res\" -type f -delete"
				    # eval "find ../Testing/SPKR1_CAFE/res/ -name \"*.res\" -type f -delete"
				    # eval "find ../Testing/SPKR1_WHITE/res/ -name \"*.res\" -type f -delete"
				    
				    eval "find ../device/ -name \"*.conf\" -type f -delete"
				    eval "find ../device/ -name \"*.phn\" -type f -delete"

				    eval "find ../MFCCs/ -name \"*.txt\" -type f -delete"
				    # eval "echo ${t} ${mnm};"
				    eval "./dtw.exe paa ${a} window ${i} banks ${j} paa_op 0 dtw_window ${n} interval_div ${d} nfft ${p} trunc ${t} mfccs 99999 test_iter 1 SIMPLE KNN knn ${k} frame_limit ${v} group_k 7 voice_k 7 SPKR1 EXPORT;"
				done
			    done
			done
		    done
		done
	    done
	done
    done
done





