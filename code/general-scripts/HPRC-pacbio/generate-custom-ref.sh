#!/bin/bash
#$ -cwd
#$ -t 1-53

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
alleledir=$prdm9dir/indiv-standard-alleles/flank10k
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
bamdir=$basedir/GRCh38-bams
customdir=$basedir/mapped-PRDM9
refdir=$customdir/refs

source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names.txt)

first=$(head -1 $basedir/gblr-calls/$sample-gblr-GRCh38.tsv | cut -f1 | cut -f1 -d/ | sed -e 's/Novel_Similar_//g' -e 's/_or_Different_Novel_Similar_.*//g')
second=$(head -1 $basedir/gblr-calls/$sample-gblr-GRCh38.tsv | cut -f1 | cut -f2 -d/ | sed -e 's/Novel_Similar_//g' -e 's/_or_Different_Novel_Similar_.*//g')

cat $alleledir/$first-flank10k.fa > $refdir/$sample-PRDM9-ref.fa
if [[ $first != $second ]]
then 
cat $alleledir/$second-flank10k.fa >> $refdir/$sample-PRDM9-ref.fa
fi
if [ $? -ne 0 ]; then fail "ref failed"; fi

echo success!
