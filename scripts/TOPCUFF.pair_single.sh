tophat=/data/Analysis/fanxiaoying/software/tophat-2.0.8b.Linux_x86_64/tophat
PATH=$PATH:/data/Analysis/fanxiaoying/software/bowtie2-2.1.0
reference=$8

outdir=$1
fq1=$2
fq2=$3
fq3=$4
bam=$5
unmap_stat=$6
insert=$7

mkdir $outdir
$tophat -p 6 -o $outdir $reference $fq1 $fq2
mv $outdir/accepted_hits.bam $outdir/Accepted_PE_hits.bam
mv $outdir/unmapped.bam $outdir/unmapped_PE.bam

$tophat -p 6 -o $outdir $reference $fq3
mv $outdir/accepted_hits.bam $outdir/Accepted_SE_hits.bam
mv $outdir/unmapped.bam $outdir/unmapped_SE.bam

perl -e 'print "@RG\tID:PE\tSM:PE\tLB:ga\tPL:Illumina\n@RG\tID:SE\tSM:hs\tLB:ga\tPL:Illumina\n"' > $outdir/rg.txt
samtools merge -rh $outdir/rg.txt -f $bam $outdir/Accepted_PE_hits.bam $outdir/Accepted_SE_hits.bam
rm $outdir/Accepted_PE_hits.bam
rm $outdir/Accepted_SE_hits.bam

samtools merge -rh $outdir/rg.txt -f Unmapped.$bam $outdir/unmapped_PE.bam $outdir/unmapped_SE.bam
rm $outdir/unmapped_PE.bam
rm $outdir/unmapped_SE.bam

samtools flagstat Unmapped.$bam >$unmap_stat 
samtools view -h -o $outdir/out.sam $bam
gtf=$9
cuff=/data/Analysis/fanxiaoying/software/cufflinks-2.1.1.Linux_x86_64/cufflinks
$cuff -o $outdir -p 4 -G $gtf $outdir/out.sam
rm $outdir/out.sam
