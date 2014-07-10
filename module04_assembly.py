from __future__ import division
import re
import subprocess
import cPickle as pickle
import time
import numpy as np
import os

class trinity(dict):
	def __init__(self,info):
		self['g'] = {}
		self['t'] = {}
		self['m0'] = info
	def run_trinity_BLAT_command(self,folder,file,reference,num,rec):
		count = 0
		filelists = []
		handle = []
		processes = []
		self['blat_cmd'] = []
		for i in range(0,num):
			filelists.append(folder+"cut_"+str(i)+'.txt')
			handle.append('f'+str(i))
			handle[i] = open(filelists[i],'w')
			processes.append("P"+str(i))
		for line in open(folder+file):
			if re.match('>',line):
				index = count%num
				count += 1
				print >>handle[index],line,
			else:
				print >>handle[index],line,
		outfile = []
		for id in range(0,num):
			i = filelists[id]
			outfile.append("out_"+i)
			handle[id].close()
			cmd = '/data/Analysis/fanxiaoying/software/blat '+reference+' '+i+' '+i+'_out.psl'
			print cmd
			self['blat_cmd'].append(cmd)
			processes[id] = subprocess.Popen(cmd,shell=True)
		time.sleep(15)
		while 1:
			count = 0
			for i in processes:
				if i.poll() is None:
					count += 1
			time.sleep(1)
			if count == 0:
				print "END"
				cmd = "cat "+folder+'*_out.psl'+'> '+folder+'/Blat.'+rec+'.psl'
				subprocess.call(cmd,shell="True")
				break
			else:
				print "RUNNING"+str(count)
	def test(self,cmd):
		subprocess.Popen(cmd,shell=True)
	def bwa_for_trinity_fa(self,fa,reads,sam,bam,statfile):
		fastafile = fa+'.fasta'
		fafile = fa+'.fa'
		bamfile = bam+'.bam'
		if not os.path.exists(fafile):
			subprocess.call(['cp',fa+'.fasta',fa+'.fa'])
		if not os.path.exists(fa+'.1.bt2'):
			cmd = '/data/Analysis/fanxiaoying/software/bowtie2-2.1.0/bowtie2-build %s %s ' %(fa+'.fa',fa)
			print "cmd"
			subprocess.call(cmd,shell=True)
		if not os.path.exists(bamfile):
			cmd1 = '/data/Analysis/fanxiaoying/software/bowtie2-2.1.0/bowtie2 -p 12 --score-min=C,-15,0 -x %s -U %s -S %s' %(fa,reads,sam)
			cmd2 = 'samtools view  -u -b -S %s | samtools sort -m 200000000 - %s' %(sam,bam)
			print "BOWTIEING"
			subprocess.call(cmd1,shell=True)
			subprocess.call(cmd2,shell=True)
			cmd3 = "samtools view -F 4 "+bamfile+" | awk '{print $3}' | perl -ne 'chomp;$hash{$_}++;END{foreach(keys %hash){print \"$_\t$hash{$_}\n\"}}' >>"+ statfile
			subprocess.call(cmd3,shell=True)
	def read_stat_file(self,samples,folder):
		self['sample'] = samples
		self['sample_info'] = {}
		for i in samples:
			self['sample_info'][i] = 0
		for i in self['g']:
			self['g'][i]['count'] = [0]*len(samples)
		for i in samples:
			sample_index = samples.index(i)
			file = open(folder+'/stat.'+i+'.txt')
			for line in file:
				array = str.split(line)
				gene = re.split('_seq',array[0])[0]
				self['g'][gene]['count'][sample_index] += int(array[1])
				self['sample_info'][i] += int(array[1])	
	def read_trinity_fa(self,file):
		fa = open(file)
		while 1:
			line = fa.readline()
			if line == '':
				break
			if line[0]=='>':
				gene = re.split(r'>|_seq',line)[1]
				trans = re.split(r'>|\s',line)[1]
				lens = int(re.split(r'len=|\s',line)[2])
				if gene not in self['g']:
					self['g'][gene] = {}
					self['g'][gene]['t'] = {}
					self['g'][gene]['longest'] = trans
					self['g'][gene]['len'] = lens
				if lens > self['g'][gene]['len']:
					self['g'][gene]['longest'] = trans
					self['g'][gene]['len'] = lens
				self['g'][gene]['t'][trans] = {}
				self['g'][gene]['t'][trans]['len'] = lens
				self['g'][gene]['t'][trans]['seq'] = ''
				self['g'][gene]['t'][trans]['genome'] = {}
				self['g'][gene]['t'][trans]['genome']['map_cond'] = 0
				self['g'][gene]['t'][trans]['genome']['maplens'] = [0]
				while len(self['g'][gene]['t'][trans]['seq']) < lens:
					seq = fa.readline()
					seq = seq.strip('\n')
					self['g'][gene]['t'][trans]['seq'] = self['g'][gene]['t'][trans]['seq']+seq
			#	print gene,trans,self['g'][gene]['t'][trans]['seq']			 
	def read_blat_file(self,filename,gen,threhold):
		file = open(filename)
		for line in file:
			array = re.split('\s*',line)
			if len(array)<10 or (not re.match('comp',array[9])):
				print array
				continue
			else:
				transc = array[9]
				gene = re.split('_seq',transc)[0]
				if gen not in self['g'][gene]['t'][transc]:
					self['g'][gene]['t'][transc][gen] = {}
					self['g'][gene]['t'][transc][gen]['map_cond'] = 0
					self['g'][gene]['t'][transc][gen]['maplens'] = []
				self['g'][gene]['t'][transc][gen]['maplens'].append(int(array[0]))
		file.close()
		for i in self['g']:
			for j in self['g'][i]['t']:
				if not gen in self['g'][i]['t'][j]:
					self['g'][i]['t'][j][gen] = {}
                                        self['g'][i]['t'][j][gen]['maplens'] = [0]
					self['g'][i]['t'][j][gen]['map_cond'] = 0
				self['g'][i]['t'][j][gen]['max'] = max(self['g'][i]['t'][j][gen]['maplens'])
				self['g'][i]['t'][j][gen]['max_90'] = len([k for k in self['g'][i]['t'][j][gen]['maplens'] if k>0.9*self['g'][i]['t'][j][gen]['max']])
		file = open(filename)
		for line in file:
			a = re.split('\s*',line)
			if len(a)<10 or (not re.match('comp',a[9])):
                                print a 
				continue
			transc = a[9]
			gene = re.split('_seq',transc)[0]
			if int(a[0]) < self['g'][gene]['t'][transc][gen]['max']:
				continue
			else:		
				self['g'][gene]['t'][transc][gen]['chr'] = a[13]
				self['g'][gene]['t'][transc][gen]['qstart'] = [int(i) for i in re.split(',',a[19]) if i != '']
				self['g'][gene]['t'][transc][gen]['tstart'] = [int(i) for i in re.split(',',a[20]) if i != '']
				self['g'][gene]['t'][transc][gen]['blockcount'] = int(a[17])
				self['g'][gene]['t'][transc][gen]['intron_50'] =  []
				self['g'][gene]['t'][transc][gen]['blocksize'] = [int(i) for i in re.split(',',a[18]) if i != '']
				self['g'][gene]['t'][transc][gen]['intronsize'] = []
				t1 = []
				t2 = []
				t3 = []
				for i in range(0,len(self['g'][gene]['t'][transc][gen]['blocksize'])):
					if self['g'][gene]['t'][transc][gen]['blocksize'][i] >= 10: 
						t1.append(self['g'][gene]['t'][transc][gen]['blocksize'][i])
						t3.append(self['g'][gene]['t'][transc][gen]['qstart'][i])
						t2.append(self['g'][gene]['t'][transc][gen]['tstart'][i])
				self['g'][gene]['t'][transc][gen]['blocksize'] = t1
				self['g'][gene]['t'][transc][gen]['tstart'] = t2
				self['g'][gene]['t'][transc][gen]['qstart'] = t3
				for i in range(1,len(self['g'][gene]['t'][transc][gen]['blocksize'])):
					self['g'][gene]['t'][transc][gen]['intronsize'].append(t2[i]-t2[i-1]-t1[i-1])
					self['g'][gene]['t'][transc][gen]['intron_50'] = [k for k in self['g'][gene]['t'][transc][gen]['intronsize'] if k>=50]
	def bedtools_merge_to_cluster(self,path,distance):
		bedfile = path+'/bed.to_merge.txt'
		sortbed = path+'/sorted_bed.to_merge.txt'
		mergebed = path +'/mergebed.txt'
		self['cluster'] = {}
		f = open(bedfile,'w')
		for gene in self['g']:
			transc =  self['g'][gene]['longest']
			info = self['g'][gene]['t'][transc]['genome']
			exon_size = len(info['blocksize'])
			print exon_size,info['tstart'],info['qstart'],info['blocksize']
			print >>f,"%s\t%s\t%s\t%s" %(info['chr'],info['tstart'][0],info['tstart'][exon_size-1]+info['blocksize'][exon_size-1],transc)
		f.close()
		cmd = "/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/bedtools sort -i %s >%s" %(bedfile,sortbed)
		print cmd
		subprocess.call(cmd,shell=True)
		cmd = "/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/bedtools merge -d %d -i %s -nms >%s" %(distance,sortbed,mergebed) 
		subprocess.call(cmd,shell=True)
		print cmd
		for line in open(mergebed):
			array = str.split(line)
			genes = re.split(r'_seq|;',array[3])
			cluster_pos = array[0]+':'+array[1]+'-'+array[2]
			if cluster_pos not in self['cluster']:
				self['cluster'][cluster_pos] = {}
				self['cluster'][cluster_pos]['count'] = [0]*len(self['sample'])
				self['cluster'][cluster_pos]['composed_genes'] = []
				self['cluster'][cluster_pos]['length'] = int(array[2])-int(array[1])
				self['cluster'][cluster_pos]['intron_50_count'] = 0
				self['cluster'][cluster_pos]['blocksize'] = []
			for i in genes:
				if not re.match('comp',i):
					continue
				if re.match('comp',i):
					self['g'][i]['cluster'] = {}
					self['g'][i]['cluster']['pos'] = array[0]+':'+array[1]+'-'+array[2]
					self['g'][i]['cluster']['range'] = int(array[2])-int(array[1])
				for j in range(0,len(self['sample'])):
					#print genes,i,self['g'][i]['count'],j,cluster_pos,self['cluster'][cluster_pos]['count']
					self['cluster'][cluster_pos]['count'][j] += self['g'][i]['count'][j]
				self['cluster'][cluster_pos]['composed_genes'].append(i)
				transc = self['g'][i]['longest']
				info = self['g'][i]['t'][transc]['genome']
				self['cluster'][cluster_pos]['intron_50_count'] += len(info['intron_50'])
				self['cluster'][cluster_pos]['blocksize'].append(info['blocksize'])
	def bedtools_merge_to_annotation(self,path,annotation_gtfs,distance,rec):
		bedfile = path+'/Trinity_Ensemble_bed.to_be_merged.txt'
		sortbed = path+'/Trinity_Ensemble_sorted_bed.to_merge.txt'
		mergebed = path +'/Trinity_Ensemble.mergebed.txt'
		f = open(bedfile,'w')
		for gene in self['g']:
                        transc =  self['g'][gene]['longest']
                        info = self['g'][gene]['t'][transc]['genome']
                        exon_size = len(info['blocksize'])
                        print >>f,"%s\t%s\t%s\t%s" %(info['chr'],info['tstart'][0],info['tstart'][exon_size-1]+info['blocksize'][exon_size-1],transc)
		for annotation_gtf in annotation_gtfs:
	                ensm = open(annotation_gtf)
			for line in ensm:
				array = re.split('\s',line)
				print >>f,"%s\t%s\t%s\t%s" %(array[0],array[1],array[2],array[3])
		f.close()
		cmd = "/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/bedtools sort -i %s >%s" %(bedfile,sortbed)
                subprocess.call(cmd,shell=True)
                cmd = "/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/bedtools merge -d %d -i %s -nms >%s" %(distance,sortbed,mergebed)
                subprocess.call(cmd,shell=True)
		for line in open(mergebed):
                        array = str.split(line)
			ENSM = []
			comp = []
			genes = re.split(r';',array[3])
			for i in genes:
				if re.findall('EN|NM|XM|NR|XR',i):
					ENSM.append(i)
				if re.findall('seq',i):
					comp.append(i)
			if comp == []:
				continue
                        for ge in comp:
				i = re.split(r'_seq',ge)[0]
                                self['g'][i][rec] = {}
                                self['g'][i][rec]['pos'] = array[0]+':'+array[1]+'-'+array[2]
                                self['g'][i][rec]['range'] = int(array[2])-int(array[1])
				self['g'][i][rec]['ensm'] = ENSM
		for i in self['cluster']:
			self['cluster'][i]['ensm'] = []
			for j in self['cluster'][i]['composed_genes']:
				for k in self['g'][j][rec]['ensm']:
					if k not in self['cluster'][i]['ensm']:
						self['cluster'][i]['ensm'].append(k)
