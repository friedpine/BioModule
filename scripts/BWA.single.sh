fq1=$2
reference=$3
outdir=$1
bam_prefix=$4
insert=$5

mkdir $outdir

bwa aln -o 1 -e 60 -i 15 -q 10 -t 8 $reference $fq1 > $outdir/$insert.1.sai
bwa samse $reference $outdir/$insert.1.sai $fq1 >$outdir/$insert.sam
samtools view  -u -b -S -t $reference.fai $outdir/$insert.sam | samtools sort -m 200000000 - $bam_prefix.SE
rm $outdir/$insert.1.sai
rm $outdir/$insert.sam
mv $bam_prefix.SE.bam $bam_prefix.bam