import subprocess,re,gzip,sys,os,pysam
bam_pe = sys.argv[1]
OP1 = gzip.open(sys.argv[2],'w')
sample = sys.argv[3]
samfile = pysam.Samfile(bam_pe,'rb')
count = 0
count_p = 0
for line in samfile:
	count += 1
	if count%10000 == 0:
		print count,count_p
	if line.is_unmapped:
		continue
	S1 = re.split('\t+',str(line))
	if line.is_read1:
		strand = 1
	else:
		strand = 2
	head = S1[0]
	seq1 = S1[9]
	OP1.write(">"+head+"::"+str(strand)+"\n"+seq1+"\n")
	count_p += 1

OP1.close()
