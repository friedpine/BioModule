#BECAREFUL: THIS USES THE PYSAM RATHER THAN READING BAMFILES DIRECTLY!!!

import re,gzip,sys,os,pysam
bam_pe = sys.argv[1]
OP1 = gzip.open(sys.argv[2],'w')
sample = sys.argv[3]
samfile = pysam.Samfile(bam_pe,'rb')
count_all = 0
count_p = 0
count_g = 0 
for line in samfile:
	count_all += 1
	if  line.is_unmapped:
		S1 = re.split('\t+',str(line))
		if S1[2] == '-1' and S1[6] == '-1':
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
			if len(seq1)>=25 and seq1.count('A')<0.6*len(seq1) and seq1.count('T')<0.6*len(seq1) and not re.findall('TATACATATGTATACATATGTATACA',seq1):
				OP1.write(">"+head+"::"+str(strand)+"\n"+seq1+"\n")
				count_g += 1
			count_p += 1
	if count_all%10000 == 0:
		print count_all,count_p,count_g
OP1.close()
