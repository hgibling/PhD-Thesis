#!/bin/bash
#$ -cwd
#$ -t 1-39
source ~/scripts/error-handling.sh
sample=NA19030
module load samtools minimap2

cd /.mounts/labs/ont/scratch/HPRC_human_data/$sample/pacbio

#cat *temp.fq > $sample-ccs.fq
#if [ $? -ne 0 ]; then fail "failed cat"; fi

#minimap2 -ax map-pb -t 20 /.mounts/labs/simpsonlab/data/references/GRCh38.fa $sample-ccs.fq > $sample-ccs.unsorted.bam
#if [ $? -ne 0 ]; then fail "failed minimap"; fi

#samtools sort $sample-ccs.unsorted.bam > $sample-ccs.bam
#if [ $? -ne 0 ]; then fail "failed samtools sort"; fi

#samtools index $sample-ccs.bam
#if [ $? -ne 0 ]; then fail "failed samtools index"; fi

#samtools view -hb $sample-ccs.bam 5:23516673-23537764 > $sample-ccs-prdm9-10k-GRCh38.bam
#if [ $? -ne 0 ]; then fail "failed samtools view"; fi

#samtools index $sample-ccs-prdm9-10k-GRCh38.bam
#if [ $? -ne 0 ]; then fail "failed samtools index 2"; fi

#samtools fastq -0 $sample-ccs-prdm9-10k.fq  $sample-ccs-prdm9-10k-GRCh38.bam
#if [ $? -ne 0 ]; then fail "failed samtools fastq"; fi

minimap2 -ax map-pb -t 20 /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/unique-PRDM9-alleles-flank10k.fa  $sample-ccs-prdm9-10k.fq > $sample-ccs-PRDM9-ref.bam
if [ $? -ne 0 ]; then fail "failed minimap 2"; fi

samtools index $sample-ccs-PRDM9-ref.bam
if [ $? -ne 0 ]; then fail "failed samtools index 3"; fi

rm $sample-ccs.fq $sample-ccs.unsorted.bam* $sample-ccs.bam*
