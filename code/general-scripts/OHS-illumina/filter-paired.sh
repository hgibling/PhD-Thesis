#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh

basedir=/.mounts/labs/simpsonlab/users/hgibling/OHS/illumina
fqdir=$basedir/PRDM9-bams/filtered-fqs

sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-sample-names.txt)

comm -3 <(cat $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_r1.fq | paste - - - - | cut -f1 | sort) <(cat $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_r2.fq | paste - - - - | cut -f1 | sort) | tr -d [:blank:] > $fqdir/$sample.filtered-unpaired.fq.names

cat $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_r1.fq | paste - - - - | sort | grep -vwf $fqdir/$sample.filtered-unpaired.fq.names | tr "\t" "\n" > $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_paired.r1.fq
cat $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_r2.fq | paste - - - - | sort | grep -vwf $fqdir/$sample.filtered-unpaired.fq.names | tr "\t" "\n" > $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_paired.r2.fq

comm -3 <(cat $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_paired.r1.fq | paste - - - - | cut -f1 | sort) <(cat $fqdir/$sample.filtered-illumina-prdm9-10k-GRCh38_paired.r2.fq | paste - - - - | cut -f1 | sort) > $fqdir/$sample-check.TEMP
if [ -s $fqdir/$sample-check.TEMP ]; then fail "$sample pairing didnt work"; fi

echo "success!"
