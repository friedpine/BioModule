import subprocess
import re
import gzip
import sys
import os
p1 = gzip.open(sys.argv[1])
OP1 = gzip.open(sys.argv[2],'wb')

total = 0
keep = 0

while 1:
	head1 = p1.readline()
	total += 1
	if head1 == '':
		INFO  = open(sys.argv[3],'w')
		print >>INFO,total/2,0,keep,0,0
		INFO.close()
		break
	seq1  = p1.readline()
	blank =	p1.readline()
	qua1  = p1.readline()
	head1 = head1.strip('\n')
        seq1  = seq1.strip('\n')
	qua1  = qua1.strip('\n')
	if not re.findall('A'*20,seq1) and not re.findall('T'*20,seq1) :
		keep += 1
		OP1.write(head1+"\n"+seq1+"\n"+"+"+"\n"+qua1+"\n")

