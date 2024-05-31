#!/bin/bash
#$ -cwd
#$ -t 1-40
source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-pacbio-links-aws.tsv | cut -f1)
file=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-pacbio-links-aws.tsv | cut -f3)
module load aws
cd /.mounts/labs/ont/scratch/HPRC_human_data
filedir=$sample/pacbio
mkdir -p $filedir
aws --no-sign-request s3 sync $file $filedir --exclude "*" --include "*.bam"
if [ $? -ne 0 ]; then fail "failed on $i"; fi
ls /.mounts/labs/ont/scratch/HPRC_human_data/$sample/pacbio/*bam > /.mounts/labs/ont/scratch/HPRC_human_data/$sample/pacbio/$sample-list.txt
