tophat=/data/Analysis/fanxiaoying/software/tophat-2.0.8b.Linux_x86_64/tophat
PATH=$PATH:/data/Analysis/fanxiaoying/software/bowtie2-2.1.0
reference=$6
outdir=$1
fq1=$2
fq2=$3
bam_PE=$4
insert=$5

mkdir $outdir
$tophat -p 6 --max-multihits 1 -o $outdir $reference $fq1 $fq2
mv $outdir/accepted_hits.bam $outdir/Accepted_PE_hits.bam
mv $outdir/unmapped.bam $outdir/unmapped_PE.bam
perl -e 'print "@RG\tID:PE\tSM:PE\tLB:ga\tPL:Illumina\n@RG\tID:SE\tSM:hs\tLB:ga\tPL:Illumina\n"' > $outdir/rg.txt
samtools merge -nrh $outdir/rg.txt -f $outdir/temp.bam $outdir/Accepted_PE_hits.bam $outdir/unmapped_PE.bam
samtools sort -n $outdir/temp.bam -@ 6 $bam_PE


#$tophat -p 6 --max-multihits 1 -o $outdir $reference $fq3
#mv $outdir/accepted_hits.bam $outdir/Accepted_SE_hits.bam
#mv $outdir/unmapped.bam $outdir/unmapped_SE.bam
#samtools merge -nrh $outdir/rg.txt -f $bam_SE $outdir/Accepted_SE_hits.bam $outdir/unmapped_SE.bam
#rm $outdir/Accepted_PE_hits.bam
#rm $outdir/Accepted_SE_hits.bam
#rm $outdir/unmapped_PE.bam
#rm $outdir/unmapped_SE.bam
