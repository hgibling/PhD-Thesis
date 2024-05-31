#!/bin/bash
#$ -cwd
#$ -t 1-39
source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-pacbio-links-aws.tsv | cut -f1)
module load samtools minimap2

cd /.mounts/labs/ont/scratch/HPRC_human_data/$sample/pacbio
pwd
echo $sample

samtools fastq -0 $sample-ccs-prdm9-10k.fq  $sample-ccs-prdm9-10k.bam
if [ $? -ne 0 ]; then fail "failed samtools fastq"; fi

minimap2 -ax map-pb -t 20 /.mounts/labs/simpsonlab/data/references/GRCh38.fa $sample-ccs-prdm9-10k.fq | samtools sort - > $sample-ccs-prdm9-10k-GRCh38.bam
if [ $? -ne 0 ]; then fail "failed minimap"; fi

samtools index $sample-ccs-prdm9-10k-GRCh38.bam
if [ $? -ne 0 ]; then fail "failed samtools index"; fi
