#!/bin/bash
#$ -cwd
source ~/scripts/error-handling.sh

cd /.mounts/labs/simpsonlab/users/hgibling/OHS/illumina/gbkc-calls

for i in *tsv
do
m=$(echo ${i##*_} | cut -f1 -d-)
s=${i%_*}
awk -v FS=',' -v OFS=',' -v m="$m" -v s="$s" '{print s, "count", m, $0}' $i >> all-methods-all-samples-gbkc-count.csv
done

wc -l all-methods-all-samples-gbkc-count.csv

