#!/bin/bash
#$ -cwd
#$ -t 1-52

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
bamdir=$basedir/GRCh38-bams
fqdir=$bamdir/fqs
customdir=$basedir/mapped-PRDM9-all

source ~/scripts/error-handling.sh
module load minimap2 samtools/1.14
sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names-used.txt)

samtools fastq -0 $fqdir/$sample-ccs-prdm9-10k-GRCh38.fq $bamdir/$sample-ccs-prdm9-10k-GRCh38.bam
if [ $? -ne 0 ]; then fail "fastq failed"; fi

minimap2 -ax map-pb $prdm9dir/all-alleles-standard-flank10k.fa $fqdir/$sample-ccs-prdm9-10k-GRCh38.fq | samtools view -hb -F 3840 | samtools sort > $customdir/$sample-ccs-prdm9-10k-GRCh38.custom-ref-all-alleles.bam
if [ $? -ne 0 ]; then fail "minimap failed"; fi
samtools index $customdir/$sample-ccs-prdm9-10k-GRCh38.custom-ref-all-alleles.bam
if [ $? -ne 0 ]; then fail "index failed"; fi

echo "success!"
