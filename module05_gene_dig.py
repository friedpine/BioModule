#!/usr/bin/python

from __future__ import division
import re
import matplotlib
matplotlib.use('Cairo')
import matplotlib.pyplot as plt
import subprocess
import cPickle as pickle
import time
import numpy as np

class genedig(dict):
	def __init__(self,cells_db):
		self['db'] = cells_db['db']
		self['sample'] = cells_db['sample']
		self['geneinfo'] = cells_db['geneinfo']
	def find_similar_genes(self,gene,normalway,rec,filename):
		self[rec] = {}
		names = []
		rvalue = []
		rules = self['geneinfo'][gene][normalway]
		for i in self['db']['lists']['ref']:
			corr = np.corrcoef(rules,self['geneinfo'][i][normalway])[1][0]
			print i,corr
			if abs(corr)>0.5:
				names.append(i)
				rvalue.append(corr)
		self[rec]['genes'] = names
		self[rec]['r'] = rvalue
		f = open(filename,'w')
		for i in range(0,len(names)):
			print >>f,names[i],rvalue[i]
			
		
