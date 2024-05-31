#!/bin/bash
#$ -cwd
#$ -t 1-51

source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
module load python samtools/1.14

gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS/

sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-samples.txt)

time $gblrdir/gap-filter.py -b $basedir/mapped-GRCh38/$sample-CCS-targeted-PRDM9.bam -v -l 84 -o $basedir/mapped-GRCh38/$sample-CCS-targeted-PRDM9-keep-reads.txt
if [ $? -ne 0 ]; then fail "gap-filter failed"; fi

samtools view -hb -N $basedir/mapped-GRCh38/$sample-CCS-targeted-PRDM9-keep-reads.txt $basedir/mapped-GRCh38/$sample-CCS-targeted-PRDM9.bam | samtools sort > $basedir/mapped-GRCh38-filtered/$sample-CCS-targeted-PRDM9-gap-filtered.bam
if [ $? -ne 0 ]; then fail "samtools view/sort failed"; fi
samtools index $basedir/mapped-GRCh38-filtered/$sample-CCS-targeted-PRDM9-gap-filtered.bam
if [ $? -ne 0 ]; then fail "samtools index failed"; fi

echo "success!"

