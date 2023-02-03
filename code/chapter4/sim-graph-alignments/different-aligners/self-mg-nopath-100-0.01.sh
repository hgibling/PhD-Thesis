#!/bin/bash
#$ -cwd
#$ -t 1-36
source ~/scripts/error-handling.sh
a=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/names.txt)
cd /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/graphs/allele-ref-graphs
module load minigraph
for i in $(seq 1 100); do
minigraph -x sr $a-nopath-chop1024.gfa /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations/READS/haploid/100/reads/$a-flank10k-100-bp-250-frag-100-x-e-0.01-i-$i.r1.fq /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations/READS/haploid/100/reads/$a-flank10k-100-bp-250-frag-100-x-e-0.01-i-$i.r2.fq | sed 's/dv:f://' | awk -v a="$a" -v i="$i" -v OFS="," '{print a, a, "nopath", "100", "0.01", i, $1, $17}' >> /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/alignments/$a-self-nopath-100-0.01-mg-deviations.csv
if [ $? -ne 0 ]; then fail "failed on $i"; fi
echo "finished $i"
done
