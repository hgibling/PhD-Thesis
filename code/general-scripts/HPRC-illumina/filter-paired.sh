#!/bin/bash
#$ -cwd
#$ -t 1-52
source ~/scripts/error-handling.sh

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-genotyping
fqdir=$basedir/fqs

sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names.txt)

comm -3 <(cat $fqdir/$sample.illumina-prdm9-10k-filtered.r1.fq | paste - - - - | cut -f1 | sort) <(cat $fqdir/$sample.illumina-prdm9-10k-filtered.r2.fq | paste - - - - | cut -f1 | sort) | tr -d [:blank:] > $fqdir/$sample.filtered-unpaired.fq.names

cat $fqdir/$sample.illumina-prdm9-10k-filtered.r1.fq | paste - - - - | sort | grep -vwf $fqdir/$sample.filtered-unpaired.fq.names | tr "\t" "\n" > $fqdir/$sample.illumina-prdm9-10k-filtered-paired.r1.fq
cat $fqdir/$sample.illumina-prdm9-10k-filtered.r2.fq | paste - - - - | sort | grep -vwf $fqdir/$sample.filtered-unpaired.fq.names | tr "\t" "\n" > $fqdir/$sample.illumina-prdm9-10k-filtered-paired.r2.fq

comm -3 <(cat $sample.illumina-prdm9-10k-filtered-paired.r1.fq | paste - - - - | cut -f1 | sort) <(cat $sample.illumina-prdm9-10k-filtered-paired.r2.fq | paste - - - - | cut -f1 | sort) > $fqdir/$sample-check.TEMP
if [ -s $fqdir/$sample-check.TEMP ]; then fail "$sample pairing didnt work"; fi

echo "success!"
