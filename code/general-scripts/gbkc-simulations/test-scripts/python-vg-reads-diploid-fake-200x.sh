#!/bin/bash
#$ -cwd

# load modules
module load python/2.7.11
module load vg
source /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/hmm-venv/bin/activate
source ~/scripts/error-handling.sh

# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling
prdm9dir=$basedir/PRDM9-Project/HMM
codedir=$prdm9dir/PRDM9-Allele-Calling/code
scratchdir=/scratch2/groups/simpsonlab/hgibling
readdir=$scratchdir/counts-haploid-python/100/reads-vg
savedir=$scratchdir/counts-diploid-python-200x/

g=$1
a1=$(echo $g | cut -d- -f1)
a2=$(echo $g | cut -d- -f2)

for i in $(seq 1 25)
do

    # combine reads
    if [ ! -s $savedir/reads-vg/$g-flank1k-100-bp-250-frag-100-x-e-0-i-$i.reads ]
    then
        cat $readdir/$a1-flank1k-100-bp-250-frag-100-x-e-0-i-$i.reads $readdir/$a2-flank1k-100-bp-250-frag-100-x-e-0-i-$i.reads > $savedir/100/reads/$g-flank1k-100-bp-250-frag-100-x-e-0-i-$i.reads

        if [ ! -s $savedir/100/reads/$g-flank1k-100-bp-250-frag-100-x-e-0-i-$i.reads ]; then fail "reads not generated for iteration $i"; fi
    fi


    echo "reads simulated for iteration $i"
done
