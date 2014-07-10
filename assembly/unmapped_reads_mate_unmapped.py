import subprocess,re,gzip,sys,os,pysam
bam_pe = sys.argv[1]
OP1 = gzip.open(sys.argv[2],'w')
sample = sys.argv[3]
samfile = pysam.Samfile(bam_pe,'rb')
count = 0
for line in samfile:
	if line.flag & 4 == 0:
		continue
	S1 = re.split('\t+',str(line))
	if S1[6] != '-1':
		continue
	if line.is_read1:
		strand = 1
	else:
		strand = 2
	head = S1[0]
	a = re.findall('#+$',S1[10]) 
	if a != []:
		seq1 = S1[9][:-len(a[0])]
	else:
		seq1 = S1[9]
	if len(seq1)>=25 and seq1.count('A')<0.6*len(seq1) and seq1.count('T')<0.6*len(seq1):
		OP1.write(">"+head+"::"+str(strand)+"\n"+seq1+"\n")
	count += 1
	if count%10000 == 0:
		print count	
OP1.close()
