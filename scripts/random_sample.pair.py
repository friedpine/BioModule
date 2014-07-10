import subprocess
import re
import gzip
import sys
import os
import numpy as np

p1 = gzip.open(sys.argv[1])
p2 = gzip.open(sys.argv[2])
OP1 = gzip.open(sys.argv[3],'wb')
OP2 = gzip.open(sys.argv[4],'wb')

while 1:
        head1 = p1.readline()
	head2 = p2.readline()
        if head1 == '':
                break
        seq1  = p1.readline()
        blank = p1.readline()
        qua1  = p1.readline()
	seq2  = p2.readline()
        blank = p2.readline()
        qua2  = p2.readline()
        if np.random.uniform(0,1)>float(sys.argv[5]):
                continue
        head1 = head1.strip('\n')
        seq1  = seq1.strip('\n')
        qua1  = qua1.strip('\n')
	head2 = head2.strip('\n')
        seq2  = seq2.strip('\n')
        qua2  = qua2.strip('\n')
	OP1.write(head1+"\n"+seq1+"\n"+"+"+"\n"+qua1+"\n")
	OP2.write(head2+"\n"+seq2+"\n"+"+"+"\n"+qua2+"\n")
