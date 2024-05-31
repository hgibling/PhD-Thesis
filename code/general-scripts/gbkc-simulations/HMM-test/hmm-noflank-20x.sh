#!/bin/bash
#$ -cwd

# load modules
module load python/2.7
source ~/scripts/error-handling.sh
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/thesispy2/bin/activate

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
codedir=$basedir/PRDM9-Allele-Calling/code
alleledir=$basedir/indiv-alleles
thesisdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/
readdir=/scratch2/groups/simpsonlab/hgibling/no-flank-20x
hmmdir=$thesisdir/PRDM9-36/HMM-test
savedir=$hmmdir/no-flank-20x/HMM

i=$SGE_TASK_ID

while read a
do
python $hmmdir/call-alleles.py -a $hmmdir/all-alleles.csv -r $readdir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i.r12.reads -e 0.0000000001 -o $savedir/$a-noflank-100-bp-250-frag-20-x-e-0-i-$i-scores
if [ $? -ne 0 ]; then fail "fail $a $i"; fi
done < $hmmdir/all-allele-names.txt

echo success!

