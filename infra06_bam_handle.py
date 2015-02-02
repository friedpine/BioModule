import pysam
import os, sys, re
import MySQLdb as mb

def get_junction_reads_counts(bamfile,chrom,left,right):
	junction_len = abs(right-left)
	count = 0
	for read in bamfile.fetch(chrom,left,right):
		if [x[1] for x in read.cigar if x[0]==3 and abs(x[1]-junction_len)<40] != []:
			count += 1
	return count

def Depth_Base_Info(samples,bamfiles,position,min_quality):
	outs = {}
	for sampleid,samfile in enumerate(bamfiles):
		print sampleid,samples[sampleid]
		sites = {}
		for read in samfile.fetch(position[0],position[1],position[2]):
			start = read.pos
			if read.mapq <= min_quality:
				continue
			cigars = [x for x in read.cigar if x[0]!=1]
			poses = [0]*len(cigars)
			for x in range(1,len(poses)):
				poses[x] = poses[x-1]+cigars[x-1][1]
			for id,seg in enumerate(cigars):
				if seg[0] != 0:
					continue
				for site in range(start+poses[id],start+poses[id]+seg[1]):
					if site+1 in sites:
						sites[site+1] += 1
					else:
						sites[site+1] = 1
		outs[samples[sampleid]] = sites
	return outs

def Junction_Info(samples,bamfiles,position,min_overlap,min_quality):
	outs = {}
	for sampleid,samfile in enumerate(bamfiles):
		print sampleid,samples[sampleid]
		juncs = {}
		for read in samfile.fetch(position[0],position[1],position[2]):
			start = read.pos
			if read.mapq <= min_quality:
				continue
			cigars = [x for x in read.cigar if x[0]!=1]
			poses = [0]*len(cigars)
			for x in range(1,len(poses)):
				poses[x] = poses[x-1]+cigars[x-1][1]
			for id,seg in enumerate(cigars):
				if seg[0] == 3:
					junc_start = start+poses[id]+1
					junc_end = start+poses[id]+seg[1]
					junc_pos = position[0]+':'+str(junc_start)+'-'+str(junc_end)
					if cigars[id-1][1]<min_overlap or cigars[id+1][1]<min_overlap:
						continue
					if junc_pos in juncs:
						juncs[junc_pos] += 1
					else:
						juncs[junc_pos] = 1
		outs[samples[sampleid]] = juncs
	return outs
