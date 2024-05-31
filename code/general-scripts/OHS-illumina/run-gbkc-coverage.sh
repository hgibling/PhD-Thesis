#!/bin/bash
#$ -cwd
#$ -t 1-50
source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/simpsonlab/users/hgibling/OHS/illumina/OHS-sample-names.txt)

prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
gbkcdir=/.mounts/labs/awadallalab/private/hgibling/gbkc/src
basedir=/.mounts/labs/simpsonlab/users/hgibling/OHS/illumina
sampledir=$basedir/PRDM9-bams/fastqs
savedir=$basedir/gbkc-calls

$gbkcdir/gbkc count -a $prdm9dir/all-alleles-standard.fa -1 $sampledir/${sample}_illumina-prdm9-10k-GRCh38_r1.fq -2 $sampledir/${sample}_illumina-prdm9-10k-GRCh38_r2.fq -d -K 131 -N 10 -f $prdm9dir/flank10k.fa -m coverage -e 0.01 -c 40 -l 151 -o $savedir/${sample}_coverage-flank.tsv
if [ $? -ne 0 ]; then fail "gbkc count failed"; fi

echo "success!"

