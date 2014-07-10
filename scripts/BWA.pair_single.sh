reference=$9
fq1=$2
fq2=$3
fq3=$4
outdir=$1
bam_prefix=$5
stat=$6
summarize=$7
insert=$8

mkdir $outdir
bwa aln -o 1 -e 60 -i 15 -q 10 -t 8 $reference $fq1 > $outdir/$insert.1.sai
bwa aln -o 1 -e 60 -i 15 -q 10 -t 8 $reference $fq2 > $outdir/$insert.2.sai
bwa sampe $reference $outdir/$insert.1.sai  $outdir/$insert.2.sai  $fq1 $fq2 >$outdir/$insert.sam
samtools view  -u -b -S -t $reference.fai $outdir/$insert.sam | samtools sort -m 200000000 - $bam_prefix.PE
rm $outdir/$insert.1.sai
rm $outdir/$insert.2.sai
rm $outdir/$insert.sam

bwa aln -o 1 -e 60 -i 15 -q 10 -t 8 $reference $fq3 > $outdir/$insert.1.sai
bwa samse $reference $outdir/$insert.1.sai $fq3 >$outdir/$insert.sam
samtools view  -u -b -S -t $reference.fai $outdir/$insert.sam | samtools sort -m 200000000 - $bam_prefix.SE
rm $outdir/$insert.1.sai
rm $outdir/$insert.sam
samtools merge -f $bam_prefix.bam $bam_prefix.PE.bam $bam_prefix.SE.bam

samtools flagstat $bam_prefix.bam >$stat
samtools view -F 4 $bam_prefix.bam | awk '{print $3}' | perl -ne 'chomp;$hash{$_}++;END{foreach(keys %hash){print "$_\t$hash{$_}\n"}}' >$summarize

rm $bam_prefix.PE.bam
rm $bam_prefix.SE.bam



