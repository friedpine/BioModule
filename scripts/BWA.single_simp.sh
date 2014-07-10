reference=$1
fq1=$2
outdir=$3
bam_prefix=$4
insert=$5

mkdir $outdir

bwa aln -o 1 -e 60 -i 15 -q 10 -t 8 $reference $fq1 > $outdir/$insert.1.sai
bwa samse $reference $outdir/$insert.1.sai $fq1 >$outdir/$insert.sam
samtools view  -u -b -S -t $reference.fai $outdir/$insert.sam | samtools sort -m 200000000 - $bam_prefix
