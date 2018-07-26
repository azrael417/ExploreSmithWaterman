#!/bin/bash
#make clean
#make -j

export OMP_NUM_THREADS=16
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
#export OMP_PROC_BIND=spread
#export OMP_PLACES=cores
#export OMP_DISPLAY_ENV=verbose

fa=/global/project/projectdirs/nstaff/cookbg/sw/100k.fa
fq=/global/project/projectdirs/nstaff/cookbg/sw/100k_illumina1.fastq

for ib in 1 2 4 6 8  `seq 10 10 100`; do
#for ib in 1 2 3; do
    # -n nbatches -b batchsize
    numactl -N 0 ./SmithWaterman -f $fa -q $fq -n 2 -b $ib 1> out.txt 2> err.txt 
    # quick sanity check
    cat out.txt |  grep cigar | sort | uniq -c
    # performance
    echo 'batch size='$ib':' $(tail -n1 out.txt)
done
