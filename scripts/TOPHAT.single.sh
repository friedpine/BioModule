tophat=/data/Analysis/fanxiaoying/software/tophat-2.0.8b.Linux_x86_64/tophat
PATH=$PATH:/data/Analysis/fanxiaoying/software/bowtie2-2.1.0

outdir=$1
fq1=$2
reference=$3
bam=$4
insert=$5
single_or_multi=$6

mkdir $outdir

mkdir $outdir
if [ $single_or_multi == 'single' ];then
	$tophat -p 10 --max-multihits 1 -o $outdir $reference $fq1
fi
if [ $single_or_multi == 'multi' ];then
	$tophat -p 10 -o $outdir $reference $fq1
fi


mv $outdir/accepted_hits.bam $bam.bam
