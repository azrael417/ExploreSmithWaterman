#!/bin/bash
#make clean
#make -j

export OMP_NUM_THREADS=16
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_DISPLAY_ENV=verbose

fa=/global/project/projectdirs/nstaff/cookbg/sw/100k.fa
fq=/global/project/projectdirs/nstaff/cookbg/sw/100k_illumina1.fastq

# -n nbatches -b batchsize
./SmithWaterman -f $fa -q $fq -n 10 -b 100 > out.txt 

# quick sanity check
cat out.txt |  grep cigar | sort | uniq -c
