#!/bin/bash
#$ -cwd

module load python/2.7

source ~/scripts/error-handling.sh
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/thesispy2/bin/activate

readdir=/scratch2/groups/simpsonlab/hgibling/test-sims-fragments/no-flank-reads
basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/HMM-test
savedir=$basedir/scores

i=$SGE_TASK_ID

while read a
do
cat $readdir/$a-noflank-100-bp-250-frag-100-x-e-0-i-$i.r1.fq $readdir/$a-noflank-100-bp-250-frag-100-x-e-0-i-$i.r2.fq | paste - - - - | cut -f2 > $readdir/$a-noflank-100-bp-250-frag-100-x-e-0-i-$i.r12.reads
python $basedir/call-alleles.py -a $basedir/all-alleles.csv -r $readdir/$a-noflank-100-bp-250-frag-100-x-e-0-i-$i.r12.reads -e 0.0000000001 -o $savedir/$a-noflank-100-bp-250-frag-100-x-e-0-i-$i-scores
if [ $? -ne 0 ]; then fail "fail $a $i"; fi
done < $basedir/all-allele-names.txt
echo success!

