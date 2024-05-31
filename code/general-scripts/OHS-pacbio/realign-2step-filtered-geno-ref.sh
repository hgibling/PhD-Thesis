#!/bin/bash
#$ -cwd
#$ -t 1-50

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS/
bamdir=$basedir/mapped-GRCh38-filtered-maxima
fqdir=$bamdir/fq2
savedir=$basedir/2step-filtered-prdm9-ref-alignments
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM

source ~/scripts/error-handling.sh
module load minimap2/2.24 samtools/1.14

sample=$(awk "NR==$SGE_TASK_ID" $savedir/samples-prdm9-ref-align.tsv | cut -f1)
first=$(awk "NR==$SGE_TASK_ID" $savedir/samples-prdm9-ref-align.tsv | cut -f4)
second=$(awk "NR==$SGE_TASK_ID" $savedir/samples-prdm9-ref-align.tsv | cut -f6)

cat $prdm9dir/indiv-standard-alleles/flank10k/$first-flank10k.fa $prdm9dir/indiv-standard-alleles/flank10k/$second-flank10k.fa > $savedir/geno-refs/$sample-geno-ref.fa

minimap2 -ax map-hifi $savedir/geno-refs/$sample-geno-ref.fa $fqdir/$sample-prdm9-spanning.fq | samtools view -hb -F 3840 | samtools sort > $savedir/$sample-prdm9-spanning-geno-ref.bam
if [ $? -ne 0 ]; then fail "minimap failed"; fi
samtools index $savedir/$sample-prdm9-spanning-geno-ref.bam
if [ $? -ne 0 ]; then fail "index failed"; fi

echo "success!"
