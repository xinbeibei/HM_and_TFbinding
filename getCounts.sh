#PBS -S /bin/bash

for histone in h3k4me2 h3k27ac h3k27me3 h3k36me3 h3k4me1 h3k4me3 h3k79me2 h3k9ac h4k20me1
do
    for rep in rep1 rep2
    do
        samtools view -c -F 4 $histone\_$rep\.bam > $histone\_$rep\_count.txt
    done
done

for rep in rep1 rep2 rep3
do
    samtools view -c -F 4 h3k9me3\_$rep\.bam > h3k9me3\_$rep\_count.txt
done
