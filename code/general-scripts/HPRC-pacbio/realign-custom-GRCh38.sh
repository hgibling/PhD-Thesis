#!/bin/bash
#$ -cwd
#$ -t 1-53
source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/sample-names.txt)
module load samtools minimap2
minimap2 -ax map-pb /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/custom-reference-genomes/refs/$sample-ref.fa /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/fastqs/$sample-ccs-prdm9-10k.fq > /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/custom-reference-genomes/$sample-realigned-custom-unsorted.bam
if [ $? -ne 0 ]; then fail "minimap2 failed"; fi
samtools sort /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/custom-reference-genomes/$sample-realigned-custom-unsorted.bam > /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/custom-reference-genomes/$sample-realigned-custom.bam
if [ $? -ne 0 ]; then fail "samtools sort failed"; fi
samtools index /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/custom-reference-genomes/$sample-realigned-custom.bam
if [ $? -ne 0 ]; then fail "samtools index failed"; fi
rm /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling/custom-reference-genomes/$sample-realigned-custom-unsorted.bam
