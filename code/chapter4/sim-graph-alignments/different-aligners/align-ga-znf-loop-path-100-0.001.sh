#!/bin/bash
#$ -cwd
#$ -t 1-36
a=$(awk "NR==$SGE_TASK_ID" /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/names.txt)
cd /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/graphs
module load GraphAligner
for i in $(seq 1 100); do
GraphAligner -g znf-loop-path-chop1024.gfa -f /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations/READS/haploid/100/reads/$a-flank10k-100-bp-250-frag-100-x-e-0.001-i-$i.r1.fq -f /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations/READS/haploid/100/reads/$a-flank10k-100-bp-250-frag-100-x-e-0.001-i-$i.r2.fq -x vg -a /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/alignments/TEMP-$a-znf-loop-path-100-0.001.gaf
sed 's/dv:f://' /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/alignments/TEMP-$a-znf-loop-path-100-0.001.gaf | awk -v a="$a" -v i="$i" -v OFS="," '{print a, "znf-loop", "path", "100", "0.001", i, $1, $14}' >> /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/alignments/$a-znf-loop-path-100-0.001-ga-deviations.csv
rm /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/alignments/TEMP-$a-znf-loop-path-100-0.001.gaf
echo "finished $i"
done
