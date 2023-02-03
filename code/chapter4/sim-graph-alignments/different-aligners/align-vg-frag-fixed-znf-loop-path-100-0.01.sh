#!/bin/bash
#$ -cwd
#$ -t 1-36
a=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/names.txt)
cd /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/graphs
module load vg
for i in $(seq 1 100); do
vg map -d znf-loop-path-chop1024 -f /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations/READS/haploid/100/reads/$a-flank10k-100-bp-250-frag-100-x-e-0.01-i-$i.r1.fq -f /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations/READS/haploid/100/reads/$a-flank10k-100-bp-250-frag-100-x-e-0.01-i-$i.r2.fq -I 5000:250:50:0:1 -U -% | sed 's/dv:f://' | awk -v a="$a" -v i="$i" -v OFS="," '{print a, "znf-loop", "path", "100", "0.01", i, $1, $16}' >> /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/alignments/$a-znf-loop-path-100-0.01-vg-frag-fixed-deviations.csv
if [ $? -ne 0 ]; then fail "failed on $i"; fi
echo "finished $i"
done
