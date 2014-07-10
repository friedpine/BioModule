tophat=/data/Analysis/fanxiaoying/software/tophat-2.0.8b.Linux_x86_64/tophat
PATH=$PATH:/data/Analysis/fanxiaoying/software/bowtie2-2.1.0
reference=/data/Analysis/fanxiaoying/database/mm10/02.GRCm38/GRCm38.no_patch_73.genome

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

samtools flagstat Unmapped.$bam >$unmap_stat 
samtools view -h -o $outdir/out.sam $bam
gtf=/data/Analysis/fanxiaoying/database/mm10/03.ensembl/Mus_musculus.GRCm38.73.gtf
cuff=/data/Analysis/fanxiaoying/software/cufflinks-2.1.1.Linux_x86_64/cufflinks
$cuff -o $outdir -p 4 -G $gtf $outdir/out.sam
rm $outdir/out.sam
