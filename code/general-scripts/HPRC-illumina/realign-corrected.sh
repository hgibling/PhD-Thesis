#!/bin/bash
#$ -cwd
#$ -t 1-52
source ~/scripts/error-handling.sh

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-genotyping
savedir=$basedir/corrected-realigned-bams
refdir=/.mounts/labs/simpsonlab/data/references

sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names.txt)
module load samtools/1.14 bwa

#bwa mem -Y -K 100000000 -t 4 $refdir/GRCh38.fa $basedir/corrected-reads/$sample.illumina-prdm9-10k-filtered-paired-corrected.r1.fa $basedir/corrected-reads/$sample.illumina-prdm9-10k-filtered-paired-corrected.r2.fa > $savedir/$sample.illumina-prdm9-10k-filtered-paired-corrected-GRCh38.unsorted.sam
#if [ $? -ne 0 ]; then fail "failed bwa"; fi
samtools sort $savedir/$sample.illumina-prdm9-10k-filtered-paired-corrected-GRCh38.unsorted.sam > $savedir/$sample.illumina-prdm9-10k-filtered-paired-corrected-GRCh38.bam
samtools index $savedir/$sample.illumina-prdm9-10k-filtered-paired-corrected-GRCh38.bam
if [ $? -ne 0 ]; then fail "failed index"; fi
rm $savedir/$sample.illumina-prdm9-10k-filtered-paired-corrected-GRCh38.unsorted.sam

echo "success!"
