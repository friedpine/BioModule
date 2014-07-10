from __future__ import division
from scipy.stats.stats import pearsonr
import re
import matplotlib.pyplot as plt
import subprocess
import cPickle as pickle
import time
import numpy as np

##BUILDING_mm10_repeats_database
def build_db():
	chr_col = 5
	left_col = 6
	right_col = 7
	rec_cols = [9,10,11,12,13,14]
	HEAD = 'T'
	infile = '/data/Analysis/fanxiaoying/database/mm10/08.repeat/rmsk.txt'
	outdb = './mm10.rmask.db'

	file = open(infile)
	db = {}
	db['chrs'] = []
	db['index'] = {}
	if HEAD == 'T':
		line = file.readline()
		a = re.split('\s+',line)
		rec_names = [a[i] for i in rec_cols]
		db['infos'] = rec_names
	for line in file:
		a = re.split('\s+',line)
		chr = a[chr_col]
		left = int(a[left_col])
		right = int(a[right_col])
		if chr not in db:
			db[chr] = {}
			db['chrs'].append(chr)
			db['index'][chr] = {}
		if left not in db[chr]:
			db[chr][left] = {}
		db[chr][left][right] = [a[i] for i in rec_cols]
		left_M = int(left/1000000)
		if left_M not in db['index'][chr]:
			db['index'][chr][left_M] = []
		db['index'][chr][left_M].append(left)
	pickle.dump(db,open(outdb,'w'))
def search_db(db,chr,left,right):
	left_M = int(left/1000000)
	for l in db['index'][chr][left_M]:
		for r in db[chr][l]:
			if (left<=r & left>=l) or (right<=r & right>=l) or (l<=right & l>=left) or (r<=right & r>=left):
				return db[chr][l][r]
build_db()
