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
		for read in samfile.pileup(position[0],position[1],position[2]):
			pos.append(column.pos)
			depth.append(column.n)
		outs.append([pos,depth])
	return outs

def Depth_Data2(bamfiles,position):
	outs = []
	for samfile in bamfiles:
		sites = {}
		for read in samfile.fetch(position[0],position[1],position[2]):
			start = read.pos
			cigars = [x for x in read.cigar if x[0]!=1]
			poses = [0]*len(cigars)
			for x in range(1,len(poses)):
				poses[x] = poses[x-1]+cigars[x-1][1]
			print start,cigars,poses
			for id,seg in enumerate(cigars):
				if seg[0] != 0:
					continue
				for site in range(start+poses[id],start+poses[id]+seg[1]):
					if site in sites:
						sites[site] += 1
					else:
						sites[site] = 1
		outs.append(sites)
	return outs

def Depth_Data2_Process_for_Plot(datas,samples,points,min_segs,whole_range,concern_ranges):
	out = []
	other_ranges = [(concern_ranges[i-1][1],concern_ranges[i][0]) for i in range(1,len(concern_ranges))]
	if whole_range[0]<concern_ranges[0][0]:
		other_ranges.insert(0,(whole_range[0],concern_ranges[0][0]))
	if concern_ranges[-1][1] < whole_range[1]:
		other_ranges.append((concern_ranges[-1][1],whole_range[1]))
	total_len = whole_range[1]-whole_range[0]
	points = min(points,total_len)
	exons_len = sum([x[1]-x[0] for x in concern_ranges])
	intron_len = total_len-exons_len
	exons_ratio = max(0.66,float(exons_len)/total_len)
	sampling_ratio_for_exons = int(exons_len/(exons_ratio*points))
	sampling_ratio_for_introns = int(intron_len/(points-exons_ratio*points))
	frame_id = []
	frame_pos = []
	frame_type = []
	for exon in concern_ranges:
		sites = range(exon[0],exon[1],min(sampling_ratio_for_exons,int((exon[1]-exon[0])/)))	
		frame_pos += sites
		frame_type += [1]*len(sites)
	for exon in other_ranges:
        sites = range(exon[0],exon[1],min(sampling_ratio_for_introns,int((exon[1]-exon[0])/)))
        frame_pos += sites
        frame_type += [0]*len(sites)
	frame_pos_sort = sorted(frame_pos)
	frame_type_sort = [frame_type[frame_pos.index(x)] for x in frame_pos_sort]
	return {frame_pos:frame_pos_sort,frame_type:frame_type_sort}

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
