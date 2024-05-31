#!/bin/bash
#$ -cwd
#$ -t 1-51
source ~/scripts/error-handling.sh
module load samtools minimap2/2.24

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS
rawdir=$basedir/raw-bams
refdir=/.mounts/labs/simpsonlab/data/references

file=$(awk "NR==$SGE_TASK_ID" $rawdir/unmapped-files.txt)
sample=$(echo $file | cut -f5 -d.)

samtools fastq -0 $rawdir/fqs/$sample-CCS-targeted-PRDM9.fq $rawdir/$file
if [ $? -ne 0 ]; then fail "fastq failed"; fi
minimap2 -ax map-hifi $refdir/GRCh38.fa $rawdir/fqs/$sample-CCS-targeted-PRDM9.fq | samtools sort > $basedir/mapped-GRCh38/$sample-CCS-targeted-PRDM9.bam
if [ $? -ne 0 ]; then fail "minimap failed"; fi
samtools index $basedir/mapped-GRCh38/$sample-CCS-targeted-PRDM9.bam
if [ $? -ne 0 ]; then fail "index failed"; fi

echo "success!"

