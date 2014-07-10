tophat=/data/Analysis/fanxiaoying/software/tophat-2.0.8b.Linux_x86_64/tophat
PATH=$PATH:/data/Analysis/fanxiaoying/software/bowtie2-2.1.0
GENONE_BAM_PREFIX=$4
infasta=$2
folder=$3
reference_genome=$1
GTF=$5

$tophat -G $GTF -p 8 -o $folder $reference_genome $infasta
samtools sort -n -@ 6 $folder'/accepted_hits.bam' $GENONE_BAM_PREFIX

#/data/Analysis/fanxiaoying/software/bowtie2-2.1.0/bowtie2 -p 8 -k 20 -x $reference_refseq -f $infasta -S $tempsam
#samtools view  -u -b -S $tempsam | samtools sort -n -@ 8 -m 1G - $REFSEQ_BAM_PREFIX
