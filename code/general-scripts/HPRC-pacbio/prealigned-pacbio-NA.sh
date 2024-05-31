#!/bin/bash
#$ -cwd
#$ -t 1-5
source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-pacbio-prealigned-links-NA-aws.tsv | cut -f1)
file=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-pacbio-prealigned-links-NA-aws.tsv | cut -f3)
module load aws samtools
cd /.mounts/labs/ont/scratch/HPRC_human_data
filedir=$sample/pacbio
#mkdir -p $filedir
#aws --no-sign-request s3 sync $file $filedir
#if [ $? -ne 0 ]; then fail "failed aws"; fi
samtools view -hb $filedir/${sample}_aligned_GRCh38_winnowmap.sorted.bam chr5:23516673-23537764 > $filedir/$sample-ccs-prdm9-10k-GRCh38.bam
if [ $? -ne 0 ]; then fail "failed samtools view"; fi
if [ -s $filedir/$sample-ccs-prdm9-10k-GRCh38.bam ]; then rm $filedir/*sorted*; samtools index $filedir/$sample-ccs-prdm9-10k-GRCh38.bam; fi
