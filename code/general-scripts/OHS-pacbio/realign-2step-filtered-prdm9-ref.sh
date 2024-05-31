#!/bin/bash
#$ -cwd
#$ -t 1-52

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS/
bamdir=$basedir/mapped-GRCh38-filtered-maxima
fqdir=$bamdir/fq2
savedir=$basedir/2step-filtered-prdm-ref-alignments
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM

source ~/scripts/error-handling.sh
module load minimap2/2.24 samtools/1.14
sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-sample-names.txt)

samtools fastq -0 $fqdir/$sample-prdm9-spanning.fq $bamdir/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.bam
if [ $? -ne 0 ]; then fail "fastq failed"; fi

minimap2 -ax map-hifi $prdm9dir/all-alleles-standard-flank10k.fa $fqdir/$sample-prdm9-spanning.fq | samtools view -hb -F 3840 | samtools sort > $savedir/$sample-prdm9-spanning-prdm9-ref.bam
if [ $? -ne 0 ]; then fail "minimap failed"; fi
samtools index $savedir/$sample-prdm9-spanning-prdm9-ref.bam
if [ $? -ne 0 ]; then fail "index failed"; fi

echo "success!"
