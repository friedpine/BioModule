pair1=$1
pair2=$2
samfile=$3
bam_prefix=$4
reference=$5
para=$6
##nature_bowtie_para_1
if [ $para == 'nature' ];then
	/data/Analysis/fanxiaoying/software/bowtie2-2.1.0/bowtie2 -p 10 --reorder --mm -M20 --score-min=C,-15,0 -x $reference -1 $pair1 -2 $pair2 -S $samfile
fi
if [ $para == 'k40' ];then
	/data/Analysis/fanxiaoying/software/bowtie2-2.1.0/bowtie2 -p 10 -k 40 --score-min=C,-15,0 -x $reference -1 $pair1 -2 $pair2 -S $samfile
fi
samtools view  -u -b -S $samfile | samtools sort -n -@ 4  - $bam_prefix
rm $samfile