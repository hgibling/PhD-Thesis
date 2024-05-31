#!/bin/bash
#$ -cwd
#$ -t 1-50

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS/
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
filterdir=$basedir/mapped-GRCh38-filtered-maxima
customdir=$basedir/mapped-PRDM9-maxima
refdir=$customdir/refs

source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-samples.txt)
module load minimap2 samtools/1.14

samtools fastq -0 $filterdir/fqs/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.fq $filterdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.bam
if [ $? -ne 0 ]; then fail "fastq failed"; fi
minimap2 -ax map-pb $refdir/$sample-PRDM9-ref.fa $filterdir/fqs/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.fq | samtools view -hb -F 3840 | samtools sort > $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.custom-ref.bam
if [ $? -ne 0 ]; then fail "minimap failed"; fi
samtools index $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.custom-ref.bam
if [ $? -ne 0 ]; then fail "index failed"; fi
echo "Success!"
