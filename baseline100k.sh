#!/bin/bash
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

fa=/global/project/projectdirs/nstaff/cookbg/sw/100k.fa
fq=/global/project/projectdirs/nstaff/cookbg/sw/100k_illumina1.fastq

numactl -N 0 ./SmithWaterman -f $fa -q $fq | tee output/baseline100k.txt
