#!/bin/bash
#$ -cwd
#$ -t 1-51

basedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/long-read-calling
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
bamdir=$basedir/GRCh38-bams
customdir=$basedir/mapped-PRDM9
refdir=$customdir/refs

source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" $basedir/sample-names.txt)
module load minimap2 samtools/1.14

first=$(head -1 $basedir/gblr-calls/$sample-gblr-GRCh38.tsv | cut -f1 | cut -f1 -d/ | sed -e 's/Novel_Similar_//g' -e 's/_or_Different_Novel_Similar_.*//g')
second=$(head -1 $basedir/gblr-calls/$sample-gblr-GRCh38.tsv | cut -f1 | cut -f2 -d/ | sed -e 's/Novel_Similar_//g' -e 's/_or_Different_Novel_Similar_.*//g')

samtools view -hb -N $basedir/gblr-calls/$sample-gblr-GRCh38.tsv.$first-reads.txt $bamdir/$sample-ccs-prdm9-10k-GRCh38.bam | samtools fastq -0 $bamdir/fqs/$sample-ccs-prdm9-10k-GRCh38.$first.fq
if [ $? -ne 0 ]; then fail "fastq failed"; fi
minimap2 -ax map-pb $refdir/$sample-PRDM9-ref.fa $bamdir/fqs/$sample-ccs-prdm9-10k-GRCh38.$first.fq | samtools view -hb -F 3840 | samtools sort > $customdir/$sample-ccs-prdm9-10k-GRCh38-$first-custom.bam
if [ $? -ne 0 ]; then fail "minimap failed"; fi
samtools index $customdir/$sample-ccs-prdm9-10k-GRCh38-$first-custom.bam
if [ $? -ne 0 ]; then fail "index failed"; fi

if [[ $first != $second ]]
then 
samtools view -hb -N $basedir/gblr-calls/$sample-gblr-GRCh38.tsv.$second-reads.txt $bamdir/$sample-ccs-prdm9-10k-GRCh38.bam | samtools fastq -0 $bamdir/fqs/$sample-ccs-prdm9-10k-GRCh38.$second.fq
if [ $? -ne 0 ]; then fail "fastq2 failed"; fi
minimap2 -ax map-pb $refdir/$sample-PRDM9-ref.fa $bamdir/fqs/$sample-ccs-prdm9-10k-GRCh38.$second.fq | samtools view -hb -F 3840 | samtools sort > $customdir/$sample-ccs-prdm9-10k-GRCh38-$second-custom.bam
if [ $? -ne 0 ]; then fail "minimap2 failed"; fi
samtools index $customdir/$sample-ccs-prdm9-10k-GRCh38-$second-custom.bam
if [ $? -ne 0 ]; then fail "index2 failed"; fi
samtools merge -o $customdir/$sample-ccs-prdm9-10k-GRCh38-custom-merged.bam $customdir/$sample-ccs-prdm9-10k-GRCh38-$first-custom.bam $customdir/$sample-ccs-prdm9-10k-GRCh38-$second-custom.bam
if [ $? -ne 0 ]; then fail "merge failed"; fi
samtools index $customdir/$sample-ccs-prdm9-10k-GRCh38-custom-merged.bam
if [ $? -ne 0 ]; then fail "index3 failed"; fi
rm $customdir/$sample-ccs-prdm9-10k-GRCh38-*-custom.bam*

else
mv $customdir/$sample-ccs-prdm9-10k-GRCh38-$first-custom.bam $customdir/$sample-ccs-prdm9-10k-GRCh38-custom-merged.bam
samtools index $customdir/$sample-ccs-prdm9-10k-GRCh38-custom-merged.bam
if [ $? -ne 0 ]; then fail "index4 failed"; fi
fi

echo success!
