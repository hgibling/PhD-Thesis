#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/simpsonlab/users/hgibling/OHS/illumina/OHS-sample-names.txt)

module load samtools

gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
basedir=/.mounts/labs/simpsonlab/users/hgibling/OHS/illumina
sampledir=$basedir/PRDM9-bams
savedir=$sampledir/fastqs

samtools fastq -1 $savedir/${sample}_illumina-prdm9-10k-GRCh38_r1.fq -2 $savedir/${sample}_illumina-prdm9-10k-GRCh38_r2.fq $sampledir/${sample}_01_LB01-01-illumina-prdm9-10k-GRCh38.bam
if [ \$? -ne 0 ]; then fail "samtools fastq failed"

echo "success!"

