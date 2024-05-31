#!/bin/bash
#$ -cwd
#$ -t 1-39
source ~/scripts/error-handling.sh

listdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject
savedir=/.mounts/labs/ont/scratch/HPRC_human_data

sample=$(awk "NR==$SGE_TASK_ID" $listdir/illumina-pacbio-prealigned-links-wget.tsv | cut -f1)
file=$(awk "NR==$SGE_TASK_ID" $listdir/illumina-pacbio-prealigned-links-wget.tsv | cut -f2)

module load aws

cd $savedir
filedir=$sample/illumina
mkdir -p $filedir
cd $filedir
wget --no-check-certificate -c --retry-connrefused $file
if [ $? -ne 0 ]; then fail "failed wget"; fi

echo "success!"
