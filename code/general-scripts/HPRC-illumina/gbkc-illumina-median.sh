#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh

listdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject
sampledir=/.mounts/labs/ont/scratch/HPRC_human_data
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
savedir=$listdir/illumina-genotyping

sample=$(awk "NR==$SGE_TASK_ID" $listdir/illumina-pacbio-prealigned-links-aws.tsv | cut -f1)

$gbkcdir/gbkc count -a $prdm9dir/all-alleles-standard.fa -1 $sampledir/$sample/illumina/$sample-illumina-GRCh38-prdm9-10k.r1.fastq -2 $sampledir/$sample/illumina/$sample-illumina-GRCh38-prdm9-10k.r2.fastq -d -K 141 -N 10 -f $prdm9dir/flank10k.fa -m median -e 0.01 -c 40 -l 151 -t 4 -o $savedir/$sample-median-flank.tsv
if [ $? -ne 0 ]; then fail "gbkc count fastq failed"; fi

echo "success!"
