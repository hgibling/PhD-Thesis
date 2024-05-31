#!/bin/bash
#$ -cwd
#$ -t 1-44
source ~/scripts/error-handling.sh

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-genotyping
savedir=/.mounts/labs/ont/scratch/HPRC_human_data
refdir=/.mounts/labs/simpsonlab/data/references

sample=$(awk "NR==$SGE_TASK_ID" $basedir/HPRC-sample-names.txt)
module load samtools/1.14 bwa

bwa mem -Y -K 100000000 -t 16 $refdir/GRCh38.fa $savedir/$sample/illumina/$sample-illumina.r1.fq.gz $savedir/$sample/illumina/$sample-illumina.r2.fq.gz > $savedir/$sample/illumina/$sample-illumina-GRCh38.unsorted.sam
if [ $? -ne 0 ]; then fail "failed bwa"; fi
samtools sort $savedir/$sample/illumina/$sample-illumina-GRCh38.unsorted.sam > $savedir/$sample/illumina/$sample-illumina-GRCh38.bam
if [ $? -ne 0 ]; then fail "failed sort"; fi
samtools index $savedir/$sample/illumina/$sample-illumina-GRCh38.bam
if [ $? -ne 0 ]; then fail "failed index"; fi
rm $savedir/$sample/illumina/$sample-illumina-GRCh38.unsorted.sam

echo "success!"
