import subprocess
import re
import gzip
import sys
import os
import numpy as np

p1 = gzip.open(sys.argv[1])
OP1 = gzip.open(sys.argv[2],'wb')

while 1:
        head1 = p1.readline()
        if head1 == '':
                break
        seq1  = p1.readline()
        blank = p1.readline()
        qua1  = p1.readline()
        if np.random.uniform(0,1)>float(sys.argv[3]):
                continue
        head1 = head1.strip('\n')
        seq1  = seq1.strip('\n')
        qua1  = qua1.strip('\n')
	OP1.write(head1+"\n"+seq1+"\n"+"+"+"\n"+qua1+"\n")
