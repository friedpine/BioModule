pair1=$1
pair2=$2
samfile=$3.sam
bam_prefix=$3
reference=$4
single_or_multi=$5
#/data/Analysis/fanxiaoying/software/bowtie2-2.1.0/bowtie2 -p 10 -k 40 --score-min=C,-15,0 -x $reference -1 $pair1 -2 $pair2 -S $samfile
if [ $single_or_multi == 'nature' ];then
	/data/Analysis/fanxiaoying/software/bowtie2-2.1.0/bowtie2 -p 10 --very-sensitive --mm -M 20 --score-min=C,-15,0 -x $reference -1 $pair1 -2 $pair2 -S $samfile
fi
if [ $single_or_multi == 'normal' ];then
    /data/Analysis/fanxiaoying/software/bowtie2-2.1.0/bowtie2 -p 10 -x $reference -1 $pair1 -2 $pair2 -S $samfile
fi

if [ $single_or_multi == 'k40' ];then
	/data/Analysis/fanxiaoying/software/bowtie2-2.1.0/bowtie2 -p 10 -k 40 --score-min=C,-15,0 -x $reference -1 $pair1 -2 $pair2 -S $samfile
fi
samtools view  -u -b -S $samfile | samtools sort -n -@ 5  - $bam_prefix
rm $samfile
