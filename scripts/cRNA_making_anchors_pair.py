import subprocess,re,gzip,sys,os
bam_pe = sys.argv[1]
OP1 = gzip.open(sys.argv[2],'w')
OP2 = gzip.open(sys.argv[3],'w')
length = int(sys.argv[4])
cmd1 = "samtools view -f 1 "+bam_pe
p1 = subprocess.Popen(cmd1,shell = True,stdout=subprocess.PIPE)
while 1:
	line1 = p1.stdout.readline()
	line2 = p1.stdout.readline()
	if line1 == '':
		break
	S1 = re.split('\s+',line1)
	S2 = re.split('\s+',line2)
	head = S1[0]
	if S1[0] != S2[0]:
		print "HEAD_different",bam_pe,head
		line3 = p1.stdout.readline()
		if line3 != '':
			S1 = re.split('\s+',line3)
		if S1[0] != S2[0]:
			line4 = p1.stdout.readline()
			if line4 != '':
				S2 = re.split('\s+',line4)
			if S1[0] != S2[0]:
				continue
	if len(S1)<10 or len(S2)<10:
		print "ONE_LINE_IS_TRANCATED",'line1',line1,'line2',line2,'END'
		continue
	if int(S1[1]) & 0x4 == 0 and int(S2[1]) & 0x4 == 0:
		continue
	if int(S1[1]) & 0x4 != 0 and len(S1[9])>=60:
		seq = S1[9]
		qua = S1[10]
		mate_info = S2[2]+"#"+S2[3]+"#"+S2[5]+'#'+'2'
		mate_seq = S2[9]
		seq1 = seq[0:length]
		seq2 = seq[len(seq)-length:len(seq)]
		qua1 = qua[0:length]
		qua2 = qua[len(seq)-length:len(seq)]
		if qua2.count('#')>10:
			jing = qua2.count('#')-5
			seq2 = seq[len(seq)-length-jing:len(seq)-jing]
			qua2 = qua[len(seq)-length-jing:len(seq)-jing]
			seq = seq[0:len(seq)-jing]
			qua = qua[0:len(seq)-jing]
		if qua1.count('#')<10 and qua2.count('#')<10:
			OP1.write("@"+head+"::"+seq+"::"+mate_seq+"\n"+seq1+"\n"+"+\n"+qua1+"\n")
			OP2.write("@"+head+"::"+seq+"::"+mate_seq+"\n"+seq2+"\n"+"+\n"+qua2+"\n")
	if int(S2[1]) & 0x4 != 0 and len(S2[9])>=60:
		seq = S2[9]
		qua = S2[10]
		mate_info = S1[2]+"#"+S1[3]+"#"+S1[5]+'#'+'1'
		mate_seq = S1[9]
		seq1 = seq[0:length]
		seq2 = seq[len(seq)-length:len(seq)]
		qua1 = qua[0:length]
		qua2 = qua[len(seq)-length:len(seq)]
		if qua2.count('#')>10:
			jing = qua2.count('#')-5
			seq2 = seq[len(seq)-length-jing:len(seq)-jing]
			qua2 = qua[len(seq)-length-jing:len(seq)-jing]
			seq = seq[0:len(seq)-jing]
			qua = qua[0:len(seq)-jing]
		if qua1.count('#')<10 and qua2.count('#')<10:
			OP1.write("@"+head+"::"+seq+"::"+mate_seq+"\n"+seq1+"\n"+"+\n"+qua1+"\n")
			OP2.write("@"+head+"::"+seq+"::"+mate_seq+"\n"+seq2+"\n"+"+\n"+qua2+"\n")
OP1.close()
OP2.close()
