#!/bin/bash
#$ -cwd
#$ -t 1-51

basedir=/.mounts/labs/simpsonlab/users/hgibling/GQ-OHS
gblrdir=/.mounts/labs/awadallalab/private/hgibling/gblr
prdm9dir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM
filterdir=$basedir/mapped-GRCh38-filtered
calldir=$basedir/gblr-calls
customdir=$basedir/mapped-PRDM9
refdir=$customdir/refs

source ~/scripts/error-handling.sh
sample=$(awk "NR==$SGE_TASK_ID" $basedir/OHS-samples.txt)
module load minimap2/2.24 samtools/1.14

first=$(head -1 $calldir/$sample-gap-filtered-gblr.tsv | cut -f1 | cut -f1 -d/ | sed -e 's/Novel_Similar_//g' -e 's/_or_Different_Novel_Similar_.*//g')
second=$(head -1 $calldir/$sample-gap-filtered-gblr.tsv | cut -f1 | cut -f2 -d/ | sed -e 's/Novel_Similar_//g' -e 's/_or_Different_Novel_Similar_.*//g')

# map all reads to all alleles
samtools fastq -0 $filterdir/$sample-CCS-targeted-PRDM9-gap-filtered-all.fq $filterdir/$sample-CCS-targeted-PRDM9-gap-filtered.bam
if [ $? -ne 0 ]; then fail "fastq failed"; fi
minimap2 -ax map-hifi $prdm9dir/all-alleles-standard-flank10k.fa $filterdir/$sample-CCS-targeted-PRDM9-gap-filtered-all.fq | samtools view -hb -F 3840 | samtools sort > $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-custom-all.bam
if [ $? -ne 0 ]; then fail "minimap failed"; fi
samtools index $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-custom-all.bam
if [ $? -ne 0 ]; then fail "index failed"; fi

# map first allele reads to first allele
samtools view -hb -N $calldir/$sample-gap-filtered-gblr.tsv.$first-reads.txt $filterdir/$sample-CCS-targeted-PRDM9-gap-filtered.bam | samtools fastq -0 $filterdir/fqs/$sample-CCS-targeted-PRDM9-gap-filtered-$first.fq
if [ $? -ne 0 ]; then fail "fastq2 failed"; fi
minimap2 -ax map-hifi $refdir/$sample-PRDM9-ref.fa $filterdir/fqs/$sample-CCS-targeted-PRDM9-gap-filtered-$first.fq | samtools view -hb -F 3840 | samtools sort > $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-$first-custom.bam
if [ $? -ne 0 ]; then fail "minimap2 failed"; fi
samtools index $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-$first-custom.bam
if [ $? -ne 0 ]; then fail "index2 failed"; fi

# map second allele reads to second allele
if [[ $first != $second ]]
then 
samtools view -hb -N $calldir/$sample-gap-filtered-gblr.tsv.$second-reads.txt $filterdir/$sample-CCS-targeted-PRDM9-gap-filtered.bam | samtools fastq -0 $filterdir/fqs/$sample-CCS-targeted-PRDM9-gap-filtered-$second.fq
if [ $? -ne 0 ]; then fail "fastq3 failed"; fi
minimap2 -ax map-hifi $refdir/$sample-PRDM9-ref.fa $filterdir/fqs/$sample-CCS-targeted-PRDM9-gap-filtered-$second.fq | samtools view -hb -F 3840 | samtools sort > $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-$second-custom.bam
if [ $? -ne 0 ]; then fail "minimap3 failed"; fi
samtools index $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-$second-custom.bam
if [ $? -ne 0 ]; then fail "index3 failed"; fi
samtools merge -o $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-custom-merged.bam $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-$first-custom.bam $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-$second-custom.bam
if [ $? -ne 0 ]; then fail "merge failed"; fi
samtools index $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-custom-merged.bam
if [ $? -ne 0 ]; then fail "index4 failed"; fi
rm $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-*-custom.bam*

# rename first allele mapped bam if hom geno
else
mv $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-$first-custom.bam $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-custom-merged.bam
samtools index $customdir/$sample-CCS-targeted-PRDM9-gap-filtered-custom-merged.bam
if [ $? -ne 0 ]; then fail "index5 failed"; fi
fi

echo success!
