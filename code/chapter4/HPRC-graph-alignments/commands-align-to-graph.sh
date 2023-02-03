##### align HPRC illumina reads to the allele-stacked graph

### setup
savedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/HPRC-graphs

cd $savedir
cp $thesisdir/PRDM9-36/graphs/hg19* $thesisdir/PRDM9-36/graphs/stacked* .
vg index -x hg19-path-chop1024.gfa.xg -g hg19-path-chop1024.gfa.gcsa hg19-path-chop1024.gfa.vg
vg index -x stacked-path-chop1024.gfa.xg -g stacked-path-chop1024.gfa.gcsa stacked-path-chop1024.gfa.vg

cd $valdir
cp geno-hprc.tsv /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/HPRC-graphs/.
cp long-read-calling/sample-info.tsv /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/HPRC-graphs/.

cd $PRDM9-bams
ln -s /.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation/PanGenomeProject/illumina-genotyping/PRDM9-10k-bams/*10k.bam* .

ls *bam | cut -f1 -d. | sort -k1,1V > ../sample-names.txt

while read s
do
frag=$(samtools stats $s.illumina-prdm9-10k.bam | grep "insert size average" | cut -f3)
std=$(samtools stats $s.illumina-prdm9-10k.bam | grep "insert size standard" | cut -f3)
echo -e "$s\t$frag\t$std" >> all-samples-frag-stats.tsv
done < ../sample-names.txt

### write alignment scripts
# code/chapter4/HPRC-graph-alignments/align.sh
# code/chapter4/HPRC-graph-alignments/align-frag.sh
# code/chapter4/HPRC-graph-alignments/align-frag-fixed.sh

### analyze
thesisdir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS
workingdir=$thesisdir/PRDM9-36/HPRC-graphs
bamdir=$workingdir/PRDM9-bams
aligndir=$workingdir/alignments

while read s
do
for g in GRCh38 graph
do
sed 's/dv:f://' $s-align-$g.gaf | awktt -v s="$s" -v g="$g" '{print s, g, "none", $1, $16}' >> all-alignments-div.tsv
sed 's/dv:f://' $s-align-$g-frag.gaf | awktt -v s="$s" -v g="$g" '{print s, g, "frag", $1, $16}' >> all-alignments-div.tsv
sed 's/dv:f://' $s-align-$g-frag-fixed.gaf | awktt -v s="$s" -v g="$g" '{print s, g, "fixed", $1, $16}' >> all-alignments-div.tsv
done
done < $workingdir/sample-names.txt

sed 's/dv:f://' 

### pull out znf reads only
savedir=/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/HPRC-graphs

cd $savedir/PRDM9-bams

while read s
do
chr=$(samtools view $s.illumina-prdm9-10k.bam | head -1 | cut -f3)
samtools view -h $s.illumina-prdm9-10k.bam $chr:23526673-23527764 | samtools fastq -1 $s.znf.reads1 -2 $s.znf.reads2
cat $s.znf.reads1 | paste - - - - | cut -f1 | sed -e 's/@//' -e 's/$/\/1/' | awktt -v s="$s" '{print s, $0}' >> all-znf-reads.tsv 
cat $s.znf.reads2 | paste - - - - | cut -f1 | sed -e 's/@//' -e 's/$/\/2/' | awktt -v s="$s" '{print s, $0}' >> all-znf-reads.tsv 
done < ../sample-names.txt

### plot
# code/chapter4/figure-4.08.R       # TODO add to repo