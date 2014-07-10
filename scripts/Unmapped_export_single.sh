bam_file=$1
out_file=$2
samtools view -f 13 $bam_file | awk '{print ">"NR;print $10}' | gzip >out_file