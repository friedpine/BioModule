tophat=/data/Analysis/fanxiaoying/software/tophat-2.0.8b.Linux_x86_64/tophat
PATH=$PATH:/data/Analysis/fanxiaoying/software/bowtie2-2.1.0
reference=$8

outdir=$1
fq1=$2
fq2=$3
reference=$4
bam=$5
insert=$6
single_or_multi=$7

mkdir $outdir
if [ $single_or_multi == 'single' ];then
	$tophat -p 10 --max-multihits 1 -o $outdir $reference $fq1 $fq2
fi
if [ $single_or_multi == 'multi' ];then
	$tophat -p 10 -o $outdir $reference $fq1 $fq2
fi

mv $outdir/accepted_hits.bam $bam.bam
