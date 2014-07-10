/data/Analysis/fanxiaoying/software/bowtie2-2.1.0/bowtie2 -p 6 -x /data/Analysis/fanxiaoying/database/hg19/00.genome/genome  -1 $1 -2 $2 -S $3
samtools view  -u -b -S $3 | samtools sort -n -@ 4  - $4
rm $3

