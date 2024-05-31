#!/bin/bash
#$ -cwd
#$ -t 1-39
source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-pacbio-links-aws.tsv | cut -f1)
module load samtools minimap2
cd /.mounts/labs/ont/scratch/HPRC_human_data/$sample/pacbio
minimap2 -ax map-pb -t 20 /.mounts/labs/simpsonlab/data/references/hs37d5.fa $sample-ccs.fq | samtools sort - > $sample-ccs-full.bam
if [ $? -ne 0 ]; then fail "failed minimap"; fi
samtools index $sample-ccs-full.bam
samtools view -hb $sample-ccs-full.bam 5:23516673-23537764 > $sample-ccs-prdm9-10k.bam
if [ $? -ne 0 ]; then fail "failed samtools view"; fi
if [ -s $sample-ccs-prdm9-10k.bam ]; then rm $sample-ccs-full.bam* $sample-ccs.fq; samtools index $sample-ccs-prdm9-10k.bam; fi
