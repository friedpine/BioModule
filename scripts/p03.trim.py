import subprocess
import re
import gzip
import sys
import os
def trim(seq,qua):
	seq = re.sub(r'(GTCGACGGCGCGCCGGATCCATA|ATATCTCGAGGGCGCGCCGGATCC|ATATGGATCCGGCGCGCCGTCGAC)',"T"*24,seq)
	try:
		capture = re.search(r'((ATATCTCG|AGGGCGCG|CCGGATCC|GTCGACGG|CGCGCCGG|GATCCATA)(.{0,30})T{10,100})',seq).group(1)
		replace = 'T'*len(capture)
		seq = re.sub(capture,replace,seq)
	except:
		a = 'A'
	if len(''.join(re.findall(r'[T]{5,}|[A]{5,}',seq)))>50:
		return 0,0,0;
	elif re.search(r'T{20,}|A{20,}',seq):
		t1=sorted(re.split(r'T{20,}|A{20,}',seq),key=len,reverse = True)[0]
		q1=qua[re.search(t1,seq).start():re.search(t1,seq).end()]
		if len(t1)>=30 and q1.count('#')<len(q1)*0.5 and (t1.count('A')+t1.count('T'))<len(t1)*0.7:
			if len(t1) == len(q1):
				return 1,t1,q1
			elif len(t1)>len(q1):
				return 1,t1[0:len(q1)],q1
			else:
				return 1,t1,q1[0:len(t1)]
		elif len(re.split(r'T{15,}|A{15,}',seq))>1:
			t2=sorted(re.split(r'T{15,}|A{15,}',seq),key=len,reverse = True)[1]
			q2=qua[re.search(t2,seq).start():re.search(t2,seq).end()]
			if len(t2)>=30 and q2.count('#')<len(q2)*0.5 and (t2.count('A')+t2.count('T'))<len(t2)*0.7:
				if len(t2) == len(q2):
					return 1,t2,q2
				elif len(t2)>len(q2):
					return 1,t2[0:len(q2)],q2
				else:
					return 1,t2,q2[0:len(t2)]
			else:
				return 0,0,0
		else:
			return 0,0,0
	else:
		return 2,seq,qua
p1 = gzip.open(sys.argv[1])
p2 = gzip.open(sys.argv[2])
OP1 = gzip.open(sys.argv[3],'wb')
OP2 = gzip.open(sys.argv[4],'wb')
OPS = gzip.open(sys.argv[5],'wb')

total = 0
keep_all = 0
keep_one_l = 0
keep_non = 0
keep_one_r = 0

while 1:
	head1 = p1.readline()
	total += 1
	if head1 == '':
		INFO  = open(sys.argv[6],'w')
		print >>INFO,total,keep_all,keep_one_l,keep_one_r,keep_non
		INFO.close()
		break
	seq1  = p1.readline()
	blank =	p1.readline()
	qua1  = p1.readline()
	head2 = p2.readline()
	seq2  = p2.readline()
	blank = p2.readline()
	qua2  = p2.readline()

	head1 = head1.strip('\n')
        seq1  = seq1.strip('\n')
	qua1  = qua1.strip('\n')
	head2 = head2.strip('\n')
	seq2  = seq2.strip('\n')
	qua2  = qua2.strip('\n')

	res1 = trim(seq1,qua1)
	res2 = trim(seq2,qua2)
	if int(res1[0])>0 and int(res2[0])>0:
		keep_all += 1	
		OP1.write(head1+"\n"+res1[1]+"\n"+"+"+"\n"+res1[2]+"\n")
		OP2.write(head2+"\n"+res2[1]+"\n"+"+"+"\n"+res2[2]+"\n")
	elif res1[0]>0:
		keep_one_l += 1
		OPS.write(head1+"\n"+res1[1]+"\n"+"+"+"\n"+res1[2]+"\n")
	elif res2[0]>0:
		keep_one_r += 1
		OPS.write(head2+"\n"+res2[1]+"\n"+"+"+"\n"+res2[2]+"\n")
	else:
		keep_non += 1
