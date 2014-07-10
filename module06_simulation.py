from __future__ import division
from scipy.stats.stats import pearsonr
import re
import os
import matplotlib.pyplot as plt
import subprocess
import cPickle as pickle
import time
import numpy as np
import subprocess

def capture_coeff_molecular_counts(n,counts,coeff,figname):
	plt.figure(figsize=(10, 8), dpi=150)
	legend_pos = [8,8,8,9,9,9]
	for i in range(n):
		count = counts[i]
		ax = plt.subplot(2,3,i+1)
		data = []
		for j in range(10000):
			captured = 0
			for k in range(count):
				if np.random.uniform(0,1)<coeff:
					captured += 1
			data.append(captured)
		ax.hist(data,label=str(count))
		plt.legend(loc=legend_pos[i])		
	plt.savefig('f.'+figname+str(coeff)+'.png')
	plt.clf()
