#!/bin/bash
#$ -cwd
#$ -t 1-2
source ~/scripts/error-handling.sh

listdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject
savedir=/.mounts/labs/ont/scratch/HPRC_human_data
refdir=/.mounts/labs/simpsonlab/data/references

sample=$(awk "NR==$SGE_TASK_ID" $listdir/illumina-pacbio-prealigned-links-aws.tsv | cut -f1)
module load bwa samtools

# match configs for other samples as described in crams
bwa mem -Y -K 100000000 -t 16 $refdir/GRCh38_full_analysis_set_plus_decoy_hla.fa $savedir/$sample/illumina/$sample-illumina-r1.fastq $savedir/$sample/illumina/$sample-illumina-r1.fastq > $savedir/$sample/illumina/$sample-illumina-GRCh38.unsorted.sam
if [ $? -ne 0 ]; then fail "failed bwa"; fi

samtools sort $savedir/$sample/illumina/$sample-illumina-GRCh38.unsorted.sam > $savedir/$sample/illumina/$sample-illumina-GRCh38.bam
samtools index $savedir/$sample/illumina/$sample-illumina-GRCh38.bam
samtools view -hb $savedir/$sample/illumina/$sample-illumina-GRCh38.bam chr5:23516673-23537764 > $savedir/$sample/illumina/$sample-illumina-GRCh38-prdm9-10k.bam
if [ $? -ne 0 ]; then fail "failed samtools view"; fi
samtools index $savedir/$sample/illumina/$sample-illumina-GRCh38-prdm9-10k.bam

echo "success!"

