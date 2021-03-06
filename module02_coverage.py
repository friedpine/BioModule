from __future__ import division
import re,os,sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pysam
import d00_sample as d00
import infra00_ranges_operate as in00


def Depth_Read_Counts(samples,bamfiles,position):
	outs = []
	for samfile in bamfiles:
		outs.append(samfile.count(position[0],position[1],position[2]))
	return outs

def DepthHash2Ranges(samples,data,min_length,min_depth):
	out = {}
	for sample in samples:
		print sample
		ranges = []
		sorted_keys = sorted(data[sample].keys())
		start = sorted_keys[0]
		end = sorted_keys[-1]
		temp = ['0']*(end-start+1)
		for x in sorted_keys:
			pos_count = 0
			total_depth = 0
			for base in range(x-9,x+1):
				if base in sorted_keys and data[sample][base]>=0.5*min_depth:
					pos_count += 1
					total_depth += data[sample][base]
			if pos_count>=5 and total_depth>=10*min_depth:
				temp[x-start] = '1'
				continue
			pos_count = 0
			total_depth = 0
			for base in range(x,x+10):
				if base in sorted_keys and data[sample][base]>=0.5*min_depth:
					pos_count += 1
					total_depth += data[sample][base]
			if pos_count>=5 and total_depth>=10*min_depth:
				temp[x-start] = '1'
		strs  = "".join(temp)
		for m in re.finditer('1{1,}',strs):
			if m.end()-m.start()>min_length:
				ranges.append([m.start()+start,m.end()+start])
		out[sample] = ranges
	return out

#Returns the mean depth of range!
def Depth_Range_Mean(samples,bamfiles,position):
	outs = []
	for sampleid,samfile in enumerate(bamfiles):
		temp = 0
		for read in samfile.fetch(position[0],position[1],position[2]):
			start = read.pos
			cigars = [x for x in read.cigar if x[0]!=1]
			poses = [0]*len(cigars)
			for x in range(1,len(poses)):
				poses[x] = poses[x-1]+cigars[x-1][1]
			for id,seg in enumerate(cigars):
				if seg[0] != 0:
					continue
				for site in range(start+poses[id],start+poses[id]+seg[1]):
					if site <= position[2] and site>=position[1]:
						temp += 1
		outs.append(temp/(position[2]-position[1]))
	return outs

def Depth_Data2(samples,bamfiles,position):
	outs = {}
	for sampleid,samfile in enumerate(bamfiles):
		sites = {}
		for read in samfile.fetch(position[0],position[1],position[2]):
			start = read.pos
			cigars = [x for x in read.cigar if x[0]!=1]
			poses = [0]*len(cigars)
			for x in range(1,len(poses)):
				poses[x] = poses[x-1]+cigars[x-1][1]
			for id,seg in enumerate(cigars):
				if seg[0] != 0:
					continue
				for site in range(start+poses[id],start+poses[id]+seg[1]):
					if site in sites:
						sites[site] += 1
					else:
						sites[site] = 1
		outs[samples[sampleid]] = sites
	return outs 

def Depth_Data2_Process_transcript(datas,samples,whole_range,introns,strand):
	out = {}
	out_ranges = in00.ranges_minus(whole_range,introns,0)
	print out_ranges
	out['id'] = range(1,sum([x[1]-x[0]+1 for x in out_ranges])+1)
	out['pos'] = []
	for out_range in out_ranges:
		out['pos'] = out['pos']+range(min(out_range),max(out_range)+1)
	if strand == '+':
		out['pos'] = sorted(out['pos'])
	if strand == '-':
		out['pos'] = sorted(out['pos'])[::-1]
	out['types'] = [1]*len(out['pos'])
	for sample in samples:
		out[sample] = [0]*len(out['pos'])
		for site in datas[sample]:
			if site in out['pos']:
				out[sample][out['pos'].index(site)] = datas[sample][site]
	return out
	
def Transc_Data_to_Ranges(datas,samples,min_len,min_depth):
	ranges = {}
	for sample in samples:
		print sample

def Look_For_Downhills(datas,samples,min_len,min_depth,min_seq):
	for sample in samples:
		d = datas[sample]
		strs  = "".join([str((sum(d[x-10:x])>=sum(d[x:x+10]) and sum(d[x-10:x])>=10)*1) for x in range(10,len(d)-10)])
		strs = '0000000000'+strs
		down_hills = []
		for m in re.finditer('1{1,}',strs):
			if m.end()-m.start()>min_len:
				down_hills.append([m.start(),m.end()])
		down_hills_sep = [down_hills[i] for i in range(0,len(down_hills)-1) if down_hills[i+1][0]>=down_hills[i][1]+min_seq]	
		down_hills_sep.append(down_hills[-1])
		for i in down_hills_sep:
			print sample,i        
	
def Depth_Data2_Process_for_Plot(datas,samples,points,min_segs,whole_range,concern_ranges):
	out = {}
	other_ranges = [(concern_ranges[i-1][1]+1,concern_ranges[i][0]-1) for i in range(1,len(concern_ranges))]
	if whole_range[0]<concern_ranges[0][0]:
		other_ranges.insert(0,(whole_range[0],concern_ranges[0][0]-1))
	if concern_ranges[-1][1] < whole_range[1]:
		other_ranges.append((concern_ranges[-1][1]+1,whole_range[1]))
	total_len = whole_range[1]-whole_range[0]
	points = min(points,total_len)
	exons_len = sum([x[1]-x[0] for x in concern_ranges])
	intron_len = total_len-exons_len
	exons_ratio = max(0.66,float(exons_len)/total_len)
	sampling_ratio_for_exons = int(exons_len/(exons_ratio*points))+1
	sampling_ratio_for_introns = int(intron_len/(points-exons_ratio*points))+1
	print other_ranges,intron_len,exons_len,exons_ratio,sampling_ratio_for_exons,sampling_ratio_for_introns
	frame_id = []
	frame_pos = []
	frame_type = []
	for exon in concern_ranges:
		sites = range(exon[0],exon[1],min(sampling_ratio_for_exons,int((exon[1]-exon[0])/min_segs)))	
		frame_pos += sites
		frame_type += [1]*len(sites)
		print 1,len(sites),exon[1],exon[0]
	for exon in other_ranges:
		sites = range(exon[0],exon[1],min(sampling_ratio_for_introns,int((exon[1]-exon[0])/min_segs)))
		frame_pos += sites
		frame_type += [0]*len(sites)
		print 0,len(sites),exon[1],exon[0]
	frame_pos_sort = sorted(frame_pos)
	frame_type_sort = [frame_type[frame_pos.index(x)] for x in frame_pos_sort]
	out['id'] = range(0,len(frame_pos_sort))
	out['pos'] = frame_pos_sort
	out['types'] = frame_type_sort
	for sample in samples:
		out[sample] = [0]*len(frame_pos_sort)
		for site in datas[sample]:
			if site in frame_pos_sort:
				out[sample][frame_pos_sort.index(site)] = datas[sample][site]
	return out


def Plot_Depth_Data(samples,datas,filename):
	plt.figure(figsize=(10, 8), dpi=150)
	n = len(samples)
	for number,sample in enumerate(samples):
		ax = plt.subplot(n,1,number+1)
		print sample,datas['id'][1:100],datas[sample][1:100]
		ax.bar(datas['id'],datas[sample],label=sample)
		ax.set_xlim([0, int(max(datas['id'])/200+1)*200])
	plt.savefig(filename)
	plt.clf()

def Samples_Bam_Handles(cursor,conn,samples,bamtype):
	bamfiles = []
	for sample in samples:
		file_path = d00.get_sample_file(cursor,sample,bamtype)
		bamfiles.append(pysam.Samfile(file_path, "rb"))
	return bamfiles