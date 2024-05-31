#!/bin/bash
#$ -cwd
#$ -t 1-50

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS/
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
filterdir=$basedir/mapped-GRCh38-filtered-maxima
customdir=$basedir/mapped-PRDM9-maxima
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM

source ~/scripts/error-handling.sh
source /.mounts/labs/simpsonlab/users/hgibling/venvs/gbkc/bin/activate
sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-samples.txt)
module load minimap2 samtools/1.14

minimap2 -ax map-pb $prdm9dir/all-alleles-standard-flank10k.fa $filterdir/fqs/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.fq | samtools view -hb -F 3840 | samtools sort > $customdir/all-alleles/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.custom-ref-all-alleles.bam
if [ $? -ne 0 ]; then fail "minimap failed"; fi
samtools index $customdir/all-alleles/$sample-CCS-targeted-PRDM9-gap-filtered-maxima.custom-ref-all-alleles.bam
if [ $? -ne 0 ]; then fail "index failed"; fi
echo "Success!"
