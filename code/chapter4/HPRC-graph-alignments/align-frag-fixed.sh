#!/bin/bash
#$ -cwd
#$ -t 1-52
source ~/scripts/error-handling.sh
source ~/.bash_profile

module load vg

thesisdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS
workingdir=$thesisdir/PRDM9-36/HPRC-graphs
bamdir=$workingdir/PRDM9-bams
aligndir=$workingdir/alignments

SAMPLE=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-frag-stats.tsv | cut -f1)
FRAG=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-frag-stats.tsv | cut -f2)
STD=$(awk "NR==$SGE_TASK_ID" $bamdir/all-samples-frag-stats.tsv | cut -f3)

vg map -% -d $workingdir/hg19-path-chop1024.gfa -b $bamdir/$SAMPLE.illumina-prdm9-10k.bam -I 5000:$FRAG:$STD:0:1 -U > $aligndir/$SAMPLE-align-GRCh38-frag-fixed.gaf
if [ $? -ne 0 ]; then fail "grch failed"; fi
vg map -% -d $workingdir/stacked-path-chop1024.gfa -b $bamdir/$SAMPLE.illumina-prdm9-10k.bam -I 5000:$FRAG:$STD:0:1 -U > $aligndir/$SAMPLE-align-graph-frag-fixed.gaf
if [ $? -ne 0 ]; then fail "graph failed"; fi
echo success!
