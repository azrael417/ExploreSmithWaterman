#!/bin/bash
#make clean
#make -j

#export OMP_NUM_THREADS=16
export OMP_PROC_BIND=spread
export OMP_PLACES=threads
#export OMP_PROC_BIND=spread
#export OMP_PLACES=cores
#export OMP_DISPLAY_ENV=verbose

fa=/global/project/projectdirs/nstaff/cookbg/sw/100k.fa
fq=/global/project/projectdirs/nstaff/cookbg/sw/100k_illumina1.fastq


for t in 1 2 4 8 12 16; do
    export OMP_NUM_THREADS=$t
    # -n nbatches -b batchsize
    numactl -N 0 ./SmithWaterman -f $fa -q $fq -n 2 -b 100  1> out.txt 2> err.txt 
    # quick sanity check
    cat out.txt |  grep cigar | sort | uniq -c
    # performance
    echo $t':' $(tail -n1 out.txt)
done
