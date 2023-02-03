#!/bin/bash
#$ -cwd
source ~/scripts/error-handling.sh

thesisdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS
savedir=$thesisdir/PRDM9-36/graphs/alignments
aligndir=/scratch2/groups/simpsonlab/hgibling/alignments-pe

a=$(awk "NR==$SGE_TASK_ID" $thesisdir/PRDM9-36/allele-znf-ends.tsv | cut -f1)
end=$(awk "NR==$SGE_TASK_ID" $thesisdir/PRDM9-36/allele-znf-ends.tsv | cut -f2)


awk -F_ -v e="$end" '$2>9900 && $3<e {print $0}' $aligndir/$a-graph-alignments.csv > $savedir/$a-graph-alignments-znf.csv
if [ $? -ne 0 ]; then fail "awk no good"; fi
echo "success!"
