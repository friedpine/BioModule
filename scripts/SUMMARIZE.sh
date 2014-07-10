bam=$1
summarize=$2
para=$3
if [ $para == 'F4' ];then
samtools view -F 4 $bam | awk '{print $3}' | perl -ne 'chomp;$hash{$_}++;END{foreach(keys %hash){print "$_\t$hash{$_}\n"}}' > $summarize
fi
