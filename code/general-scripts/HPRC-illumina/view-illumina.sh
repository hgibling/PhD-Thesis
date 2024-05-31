#!/bin/bash
#$ -cwd
#$ -t 3-39
source ~/scripts/error-handling.sh

listdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject
savedir=/.mounts/labs/ont/scratch/HPRC_human_data
refdir=/.mounts/labs/simpsonlab/data/references

sample=$(awk "NR==$SGE_TASK_ID" $listdir/illumina-pacbio-prealigned-links-aws.tsv | cut -f1)
module load samtools

samtools index $savedir/$sample/illumina/$sample.final.cram
samtools view -hb $savedir/$sample/illumina/$sample.final.cram chr5:23516673-23537764 > $savedir/$sample/illumina/$sample-illumina-GRCh38-prdm9-10k.bam
if [ $? -ne 0 ]; then fail "failed samtools view"; fi
samtools index $savedir/$sample/illumina/$sample-illumina-GRCh38-prdm9-10k.bam
if [ $? -ne 0 ]; then fail "failed samtools index"; fi

echo "success!"
