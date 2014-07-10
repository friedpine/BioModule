import subprocess,re,gzip,sys,os,pysam
bam_pe = sys.argv[1]
OP1 = gzip.open(sys.argv[2],'w')
OP2 = gzip.open(sys.argv[3],'w')
length = int(sys.argv[4])
samfile = pysam.Samfile(bam_pe,'rb')
count = 0
line1 = ['','']
for line in samfile:
	line2 = re.split('\t+',str(line))
	if line2[0] == line1[0]:
		S1 = line1
		S2 = line2
		head = S1[0]
		a = re.findall('#+$',S1[10]) 
		if a != []:
			seq1 = S1[9][:-len(a[0])]
			qua1 = S1[10][:-len(a[0])]
		else:
			seq1 = S1[9]
			qua1 = S1[10]
		a = re.findall('#+$',S2[10]) 
		if a != []:
			seq2 = S2[9][:-len(a[0])]
			qua2 = S2[10][:-len(a[0])]
		else:
			seq2 = S2[9]
			qua2 = S2[10]
		if len(seq1)>50 and int(S1[1]) & 0x4 != 0 and seq1.count('A')<0.6*len(seq1) and seq1.count('T')<0.6*len(seq1):
			OP1.write("@"+head+"_1::"+seq1+"::"+seq2+"\n"+seq1[0:length]+"\n"+"+\n"+qua1[0:length]+"\n")
			OP2.write("@"+head+"_1::"+seq1+"::"+seq2+"\n"+seq1[-length:]+"\n"+"+\n"+qua1[-length:]+"\n")
		if len(seq2)>50 and int(S2[1]) & 0x4 != 0 and seq2.count('T')<0.6*len(seq2) and seq2.count('T')<0.6*len(seq2):
			OP1.write("@"+head+"_2::"+seq2+"::"+seq1+"\n"+seq2[0:length]+"\n"+"+\n"+qua2[0:length]+"\n")
			OP2.write("@"+head+"_2::"+seq2+"::"+seq1+"\n"+seq2[-length:]+"\n"+"+\n"+qua2[-length:]+"\n")	
	line1 = line2
	count += 1
	if count%10000 == 0:
		print count	
OP1.close()
OP2.close()
