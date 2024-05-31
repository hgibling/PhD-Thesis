#!/bin/bash
#$ -cwd
#$ -t 1-3
source ~/scripts/error-handling.sh

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/TRIO
savedir=$basedir/redo-corrected
refdir=/.mounts/labs/simpsonlab/data/references

sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names.txt)
module load samtools/1.14 bwa

bwa mem -Y -K 100000000 -t 4 $refdir/GRCh38.fa $savedir/$sample-300x-r1.paired $savedir/$sample-300x-r2.paired > $savedir/$sample-300x-GIAB-PRDM9-corrected.unsorted.sam
if [ $? -ne 0 ]; then fail "failed bwa"; fi
samtools sort $savedir/$sample-300x-GIAB-PRDM9-corrected.unsorted.sam > $savedir/$sample-300x-GIAB-PRDM9-corrected.sam
samtools index $savedir/$sample-300x-GIAB-PRDM9-corrected.sam
if [ $? -ne 0 ]; then fail "failed index"; fi
rm $savedir/$sample-300x*unsorted.sam

echo "success!"
