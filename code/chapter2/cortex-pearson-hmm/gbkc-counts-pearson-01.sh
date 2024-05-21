#!/bin/bash
#$ -cwd

# command use: ./gbkc-counts-scoring.sh $1 allele list $2 coverage $3 read length $4 error rate $5 fragment length $6 max kmer size $7 method $8 iteration number


source ~/scripts/error-handling.sh


# define directories
basedir=/.mounts/labs/awadallalab/private/hgibling
thesisdir=$basedir/PRDM9-Project/THESIS

savedir=$thesisdir/PRDM9-36/pearson-test/scores
simdir=$basedir/PRDM9-Project/SIMULATIONS/haploid
gbkcdir=$basedir/gbkc/src


# make directories if they don't already exist

r=100
e=0.01
c=100
f=250
k=99
m=coverage

i=$SGE_TASK_ID


# iterate over list of alleles
while read a
do
    for c in 20
    do
        # calculate scores
        if [ ! -s $savedir/$a-flank10k-$r-bp-$f-frag-$c-x-e-$e-i-$i-$m-m-scores-full.csv ]
        then
            $gbkcdir/gbkc count -a $thesisdir/PRDM9-36/all-alleles.fa -1 $simdir/$c/reads/$a-flank10k-$r-bp-$f-frag-$c-x-e-$e-i-$i.r1.fq -2 $simdir/$c/reads/$a-flank10k-$r-bp-$f-frag-$c-x-e-$e-i-$i.r2.fq -f $thesisdir/flank10k.fa -K $k -l $r -e $e -c $c -m $m -o $savedir/$a-flank10k-$r-bp-$f-frag-$c-x-e-$e-i-$i-$m-m-scores.csv
            if [ $? -ne 0 ]; then fail "gbkc failed at allele $a"; fi
            awk -v a="$a" '{print a "," $0}' $savedir/$a-flank10k-$r-bp-$f-frag-$c-x-e-$e-i-$i-$m-m-scores.csv > $savedir/$a-flank10k-$r-bp-$f-frag-$c-x-e-$e-i-$i-$m-m-scores-full.csv
            if [ ! -s $savedir/$a-flank10k-$r-bp-$f-frag-$c-x-e-$e-i-$i-$m-m-scores-full.csv ]; then fail "No scores for allele $a"; else rm $savedir/$a-flank10k-$r-bp-$f-frag-$c-x-e-$e-i-$i-$m-m-scores.csv; fi
        fi
    done
	echo "Scores calculated for allele $a"
done < $thesisdir/PRDM9-36/all-allele-names.txt

echo "success!"
