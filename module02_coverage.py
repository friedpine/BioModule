from __future__ import division
import re
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pysam
import d00_sample as d00


def Depth_Data(bamfiles,position):
	outs = []
	for samfile in bamfiles:
		pos = []
		depth = []
		for read in samfile.pileup(posit ion[0],position[1],position[2]):
			pos.append(column.pos)
			depth.append(column.n)
		outs.append([pos,depth])
	return outs

def Depth_Data2(bamfiles,position):
	outs = []
	for samfile in bamfiles:
		reads = []
		for read in samfile.fetch(posit ion[0],position[1],position[2]):
			reads.append(read)
		outs.append(reads)
	return outs
	

def Plot_Depth_Data(samples,depths,filename):
	plt.figure(figsize=(10, 8), dpi=150)
	n = len(samples)
	for i in range(n):
		ax = plt.subplot(n,1,i+1)
		ax.bar(range(len(depths[1])),depths[i],label=samples[i])
		leg = plt.legend(2)
		leg.draw_frame(False)
	plt.savefig(filename)
	plt.clf()

def Samples_Bam_Handles(cursor,conn,samples,bamtype):
	bamfiles = []
	for sample in samples:
		file_path = d00.get_sample_file(cursor,sample,bamtype)
		bamfiles.append(pysam.Samfile(file_path, "rb"))
	return bamfiles

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
