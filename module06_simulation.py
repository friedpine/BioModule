from __future__ import division
import re
import os
import subprocess
import time
import numpy as np
import d00_sample as d00

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

def fq_samplings(cursor,conn,dbname,samples,readcounts,fq1_types,fq2_types,outdir,file):
	f = open(file,"wb")
	for id,sample in enumerate(samples):
		fq1 = ""
		fq2 = ""
		try:
			for i,j in enumerate(fq1_types):
				fq1 = d00.get_path2(cursor,sample,j)
				fq2 = d00.get_path2(cursor,sample,fq2_types[i])
		except:
			tmp = 1
		if fq1 == "" or fq2 == "":
			print "NO_FQ_FILES!!"
		else:
			tmp_fq1 = outdir+"/"+sample+"_1.fq"
			tmp_fq2 = outdir+"/"+sample+"_2.fq"
			print >>f, "less "+fq1+" | head -n "+str(readcounts*4)+" >"+tmp_fq1
			print >>f, "less "+fq2+" | head -n "+str(readcounts*4)+" >"+tmp_fq2
			d00.insert_sample_file(cursor,conn,sample,fq1+"_"+str(readcounts),tmp_fq1)
			d00.insert_sample_file(cursor,conn,sample,fq2+"_"+str(readcounts),tmp_fq2)
			
			
	
