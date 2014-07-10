tophat=/data/Analysis/fanxiaoying/software/tophat-2.0.8b.Linux_x86_64/tophat
PATH=$PATH:/data/Analysis/fanxiaoying/software/bowtie2-2.1.0
reference=/data/Analysis/fanxiaoying/database/mm10/02.GRCm38/GRCm38.no_patch_73.genome

outdir=$1
fq1=$2
fq2=$3
fq3=$4
bam_PE=$5
bam_SE=$6
insert=$7
threads=$8

mkdir $outdir
#$tophat -p 6 --max-multihits 1 -o $outdir $reference $fq1 $fq2
#mv $outdir/accepted_hits.bam $outdir/Accepted_PE_hits.bam
#mv $outdir/unmapped.bam $outdir/unmapped_PE.bam
perl -e 'print "@RG\tID:PE\tSM:PE\tLB:ga\tPL:Illumina\n@RG\tID:SE\tSM:hs\tLB:ga\tPL:Illumina\n"' > $outdir/rg.txt
#samtools merge -nrh $outdir/rg.txt -f $bam_PE $outdir/Accepted_PE_hits.bam $outdir/unmapped_PE.bam

$tophat -p $threads --max-multihits 1 -o $outdir $reference $fq3
mv $outdir/accepted_hits.bam $outdir/Accepted_SE_hits.bam
mv $outdir/unmapped.bam $outdir/unmapped_SE.bam
samtools merge -nrh $outdir/rg.txt -f $bam_SE $outdir/Accepted_SE_hits.bam $outdir/unmapped_SE.bam
#rm $outdir/Accepted_PE_hits.bam
#rm $outdir/Accepted_SE_hits.bam
#rm $outdir/unmapped_PE.bam
#rm $outdir/unmapped_SE.bam
