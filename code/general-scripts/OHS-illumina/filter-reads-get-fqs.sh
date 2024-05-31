#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh
module load samtools

gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
basedir=/.mounts/labs/simpsonlab/users/hgibling/OHS/illumina
sampledir=$basedir/PRDM9-bams
savedir=$sampledir/filtered-fqs

sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-sample-names.txt)

samtools view -hb -F 3844 $sampledir/${sample}_01_LB01-01-illumina-prdm9-10k-GRCh38.bam > $sampledir/${sample}.filtered-illumina-prdm9-10k-GRCh38.bam
if [ $? -ne 0 ]; then fail "view failed"; fi
samtools index $sampledir/${sample}.filtered-illumina-prdm9-10k-GRCh38.bam
if [ $? -ne 0 ]; then fail "index failed"; fi
samtools fastq -1 $savedir/${sample}.filtered-illumina-prdm9-10k-GRCh38_r1.fq -2 $savedir/${sample}.filtered-illumina-prdm9-10k-GRCh38_r2.fq $sampledir/${sample}.filtered-illumina-prdm9-10k-GRCh38.bam
if [ $? -ne 0 ]; then fail "fastq failed"; fi

echo "success!"

