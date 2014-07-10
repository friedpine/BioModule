from __future__ import division
from scipy.stats.stats import pearsonr
import re
import os
import matplotlib
matplotlib.use('Cairo')
import matplotlib.pyplot as plt
import subprocess
import cPickle as pickle
import time
import numpy as np
import subprocess


def get_coverage(bamfile,ID,length):
	cmd = "samtools view -F 4 "+bamfile+" "+ID
	p1 = subprocess.Popen(cmd,shell = True,stdout=subprocess.PIPE)
	depth = [0]*length
	for line in p1.stdout:
		array = str.split(line)
		pos = int(array[3])
		info = re.findall(r'(\d+)(\D+)',array[5])
		for i in range(0,len(re.findall(r'(\d+)(\D+)',array[5]))):
			A = int(info[i][0])
			B = info[i][1]
			if B == 'M':
				for j in range(pos-1,pos+A-1):
					try:
						depth[j] += 1
					except:
						HEHE = 'Hehe'
				pos += A
			elif B == "I":
				pos += 0
			else:
				pos += A
	return depth
class gene_coverage_depth(dict):
	def __init__(self):
                self['coverage'] = {}
	def get_gene_coverage(self,geneclass,gene,db1,db2,samples):
		transcs = []
		all_exon_pos = []
		self['coverage']['merge_exons'] = 'AAAA'
		merge_exon_pos = []
		if gene not in db1[geneclass]['g']:
			print "GENE NOT EXISTS IN DB",gene
			return 0
		for i in db1[geneclass]['g'][gene]['transcript']:
			transcs.append(i)
			for j in range(0,int(len(db2['transc_info'][i]['exons_pos'])/2)):
				exon_pos = [db2['transc_info'][i]['exons_pos'][2*j],db2['transc_info'][i]['exons_pos'][2*j+1]]
				if exon_pos not in all_exon_pos:
					all_exon_pos.append(exon_pos)
		for i in range(0,len(all_exon_pos)):
			for j in range(i,len(all_exon_pos)):
				if (all_exon_pos[j][0]>= all_exon_pos[i][0] and all_exon_pos[j][0]<= all_exon_pos[i][1]) or (all_exon_pos[j][1] >= all_exon_pos[i][0] and all_exon_pos[j][1] <= all_exon_pos[i][1]) or (all_exon_pos[i][0]>= all_exon_pos[j][0] and all_exon_pos[i][0]<= all_exon_pos[j][1]):
					print all_exon_pos[i],all_exon_pos[j],
					m_left = min(all_exon_pos[j][0],all_exon_pos[i][0])
					m_right = max(all_exon_pos[j][1],all_exon_pos[i][1])
					all_exon_pos[j] = [m_left,m_right]
					print all_exon_pos[j]
					all_exon_pos[i]	= [0,0]
		self['coverage']['merge_exons'] = sorted([i for i in all_exon_pos if i!=[0,0]])
		print self['coverage']['merge_exons'],all_exon_pos
		for i in all_exon_pos:
			print i
class coverage(dict):
	def initialize(self,database):
		self['db'] = database
	def build_coverage(self,cat,counts_file,bamfile):
		count = 0
		count_well = 0
		self[cat] = {}
		self[cat]['well_expressed_list'] = []
		file = open(counts_file)
		for line in file:
			count += 1
			array = str.split(line)
			trans = array[0]
			try:		
				gene = self['db']['ref']['t'][array[0]]
			except:
				continue
			print trans,gene,count,
			length = self['db']['ref']['g'][gene]['transcript'][array[0]]
			if 100*int(array[1])/length>=10 and self['db']['ref']['g'][gene]['counts']==1:
				count_well += 1
				print count_well,
				self[cat]['well_expressed_list'].append(trans)
				self[cat][trans] = {}
				self[cat][trans]['depth'] = get_coverage(bamfile,trans,length)
				self[cat][trans]['raw'] = [0]*20
				for i in range(0,20):
					total = sum(self[cat][trans]['depth'])
					self[cat][trans]['raw'][i] = max(sum(self[cat][trans]['depth'][int(length/20)*i:int(length/20)*(i+1)])/total,0.0005)
	def normalize_by_bulk(self,bulk,sample):
		self[sample]['normal_list'] = []
		for i in self[sample]['well_expressed_list']:
			if i in self[bulk]:
				print "##########",i,
				self[sample]['normal_list'].append(i)
				self[sample][i]['normal'] = [0]*20
				for j in range(0,20):
					if self[bulk][i]['raw'][j]<0.004:
						self[sample][i]['normal'][j] = 'na'
					else:
						self[sample][i]['normal'][j] = self[sample][i]['raw'][j]/self[bulk][i]['raw'][j]
				print self[sample][i]['normal']
	def bin_by_length(self,sample,geneslist,normal_way):
		file = "length."+sample+".txt"
		f = open(file,'w')
		self[sample]['bin'] = [0,1,2,3,4,5]
		self[sample]['bin_std'] = [0,1,2,3,4,5]
		temp = [0,1,2,3,4,5]
		for i in range(0,6):
			self[sample]['bin'][i] = [0]*20
			self[sample]['bin_std'][i] = [0]*20
			temp[i] = [0]*20
			for j in range(0,20):
                                temp[i][j] = []
		for i in self[sample][geneslist]:
			kbin = min(int(len(self[sample][i]['depth'])/1000),5)
			for j in range(0,20):
				if self[sample][i][normal_way][j] != 'na':
					print kbin,i,j
					temp[kbin][j].append(self[sample][i][normal_way][j])
		for i in range(0,6):
			for j in range(0,20):
				self[sample]['bin'][i][j] = np.average(temp[i][j])
				self[sample]['bin_std'][i][j] = np.std(temp[i][j],ddof=1)
		print >>f,self[sample]['bin']
		print >>f,self[sample]['bin_std']
	def plot_all_six_lines(self,samples,labels,legends,colors,record):
		plt.figure(figsize=(10, 8), dpi=150)
		legend_pos = [8,8,8,9,9,9]
                for i in range(6):
			ax = plt.subplot(2,3,i+1)
			for j in range(len(samples)):
				print self[samples[j]]['bin'][i],legends[j],colors[j]
				ax.plot(self[samples[j]]['bin'][i],linewidth=2,linestyle="-",label=legends[j],c = colors[j])
				ax.set_xticks([0,5,10,15,20])
				ax.set_xticklabels([0,25,50,75,100])
			leg = plt.legend(loc=legend_pos[i])
			leg.draw_frame(False)
			plt.title(labels[i])
		plt.savefig('f.'+record+'.png')
                plt.clf()
	def plot_six_lines_with_reading_files(self,samples,colors,trans,labels,range_name,record):
		length_raw = {}
		length = {}
		count = 0
		for i in samples:
			length[labels[count]] = {}
			for s in i:
				length_raw[s] = {}
				f = open('./length/length.'+s+'.txt')
				line = f.readline()
				line_split = re.split('\[|\]|,',line)
				temp = []
				for j in line_split:
					if '0' in j:
						temp.append(float(j))
				length_raw[s][1] = temp[0:20]
				length_raw[s][2] = temp[20:40]
				length_raw[s][3] = temp[40:60]
				length_raw[s][4] = temp[60:80]
				length_raw[s][5] = temp[80:100]
				length_raw[s][6] = temp[100:120]
			for m in range(1,7):
				length[labels[count]][m] = {}
	                        length[labels[count]][m]['mean'] = []
	                        length[labels[count]][m]['std'] = []
				for n in range(0,20):
					temp = []
					for p in i:
						print i,p,m,n,length[p][m]
						temp.append(length_raw[p][m][n])
					length[labels[count]][m]['mean'].append(np.mean(temp))
					length[labels[count]][m]['std'].append(np.std(temp,ddof=1))
			count += 1
		plt.figure(figsize=(10, 8), dpi=200)
		for i in range(5):
			ax = plt.subplot(2,3,i+1)
			for j in range(len(labels)):
				print range(0,20)
				print length[labels[j]][i+1]['mean']
				print length[labels[j]][i+1]['std']
				ax.errorbar(range(0,20),length[labels[j]][i+1]['mean'],yerr=length[labels[j]][i+1]['std'],linewidth=2,linestyle="-",c = colors[j],alpha = trans[j])
				ax.set_xticks([0,5,10,15,20])
				ax.set_xticklabels([0,25,50,75,100])
			plt.title(range_name[i])
		i = 5
		for j in range(len(labels)):
			ax = plt.subplot(2,3,i+1)
			ax.errorbar(range(0,20),length[labels[j]][i+1]['mean'],yerr=length[labels[j]][i+1]['std'],linewidth=2,linestyle="-",c = colors[j],alpha = trans[j],label=labels[j])
	       	        ax.set_xticks([0,5,10,15,20])
	                ax.set_xticklabels([0,25,50,75,100])
		plt.ylim(0,0.30)
		plt.title(range_name[i])
		plt.legend(loc=2)	
		plt.savefig('f.png.'+record+'.png',dpi=200)
		plt.savefig('f.eps.'+record+'.eps')
                plt.clf()
	def GC_content_build_bulk(self,n,GC_file,bulk_sample):
                self['GC'] = {}
		self['GC']['gene'] = {}
                for i in range(n):
                        self['GC'][i] = {}
			self['GC'][i]['bulk'] = []
                bulk = open(bulk_sample)
                GC_file = open(GC_file)
                for line in GC_file:
                        b = re.split('\s',line)
                        #GC_100 = int(float(b[3])*n)
                        #self['GC']['gene'][b[0]] = GC_100
			self['GC']['gene'][b[0]] = min(int(int(b[1])/500),9)
		for line in bulk:
			b = re.split('\s',line)
			i = self['GC']['gene'][b[0]]
			self['GC'][i]['bulk'].append(b[0])
	def GC_content_sample(self,n,samples,labels):
		count = 0
		for sample in samples:
			for i in range(n):
				self['GC'][i][labels[count]] = []
			self['GC'][labels[count]] = []
			f = open(sample)
			for line in f:
				b = re.split('\s',line)
				GC_c = self['GC']['gene'][b[0]] 
				if b[0] in self['GC'][GC_c]['bulk']:
					self['GC'][GC_c][labels[count]].append(b[0])
			count += 1	
