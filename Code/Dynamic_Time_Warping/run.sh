#!/bin/bash
for i in 64 128 # ...
do
    z=$(($i-(($i/10))))
    r=$(($i/5)) # +(($i/4))))
    for((j=32;j<=z;j=j+$((i/10))))
    do
	for ((n=100;n<=1000;n=n+200)) # n in 200  #in 0 $((1/10)) # 16, 32, etc.
	do
	    for p in 256 # $i and 32
	    do
		for d in 2
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

		    eval "./dtw.exe paa 2 window ${i} banks ${j} paa_op 0 dtw_window ${n} interval_div ${d} nfft ${p} mfccs 999 FEMALE SIMPLE;"
		done
	    done
	done
    done
done
