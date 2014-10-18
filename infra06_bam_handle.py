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



