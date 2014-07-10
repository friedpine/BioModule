##lists: all kinds of subgroup genelists
##lists_ids: the index of the subgroup genelists
##db['array_gene']: the ORDER all/whole genesets!! IN ORDER!!
##db['genes']: the whole gene BASIC information!!!
##db['ens']: store the whole ensemble information!!!
##db['ref']: store the whole refseq information!!!
from __future__ import division
from scipy.stats.stats import pearsonr
import re,sys,os,subprocess,time
import matplotlib.pyplot as plt
import cPickle as pickle
import numpy as np


db = {}
db['lists'] = {}
db['lists_ids'] = {}
db['ref'] = {}
db['ens'] = {}
db['ref']['g'] = {}
db['ref']['t'] = {}
db['ens']['g'] = {}
db['ens']['t'] = {}

file = open('/data/Analysis/fanxiaoying/database/mm9/09.python_db_make/info.txt')
for line in file:
	array = str.split(line)
	if array[1] not in db['ref']['g']:
		db['ref']['g'][array[1]] = {}
		db['ref']['g'][array[1]]['transcript'] = {}
		db['ref']['g'][array[1]]['counts'] = 0
	db['ref']['g'][array[1]]['transcript'][array[0]] = int(array[2])
	db['ref']['g'][array[1]]['counts'] += 1
	db['ref']['t'][array[0]] = array[1]
for i in db['ref']['g']:
	db['ref']['g'][i]['length'] = 100
	for j in db['ref']['g'][i]['transcript']:	
		if db['ref']['g'][i]['transcript'][j] > db['ref']['g'][i]['length']:
			db['ref']['g'][i]['length'] = db['ref']['g'][i]['transcript'][j]

###READ_THE_ENSEMBL_GENES
file = open('/data/Analysis/fanxiaoying/database/mm9/09.python_db_make/ensmbl.info.txt')
for line in file:
	array = str.split(line)
	if array[1] not in db['ens']['g']:
		db['ens']['g'][array[1]] = {}
		db['ens']['g'][array[1]]['class'] = array[4]
		db['ens']['g'][array[1]]['ENSMUSG'] = array[4]
		db['ens']['g'][array[1]]['transcript'] = {}
		db['ens']['g'][array[1]]['counts'] = 0
	db['ens']['g'][array[1]]['transcript'][array[0]] = int(array[2])
	db['ens']['g'][array[1]]['counts'] += 1
	db['ens']['t'][array[0]] = array[1]
for i in db['ens']['g']:
	db['ens']['g'][i]['length'] = 100
	for j in db['ens']['g'][i]['transcript']:
		if db['ens']['g'][i]['transcript'][j] > db['ens']['g'][i]['length']:
			db['ens']['g'][i]['length'] = db['ens']['g'][i]['transcript'][j]
###build_genes:
db['genes'] = {}
db['array_gene'] = []

###GENELIST:REFSEQ_AND_BUILD_GENES_INFO
db['lists']['ref'] = []
db['lists']['ERCC'] = []
db['lists_ids']['ref'] = []
db['lists_ids']['ERCC'] = []
for i in db['ref']['g']:
	db['genes'][i] = {}
	db['array_gene'].append(i)
	if re.findall('ERCC',i):
		db['lists']['ERCC'].append(i)
		db['lists_ids']['ERCC'].append(len(db['array_gene'])-1)
		db['genes'][i]['reference'] = 'ref'
	else:
		db['lists']['ref'].append(i)
		db['lists_ids']['ref'].append(len(db['array_gene'])-1)
		db['genes'][i]['reference'] = 'ref'
	db['genes'][i]['length'] = db['ref']['g'][i]['length']

###GENELIST:ENSEMBL_AND_BUILD_GENES_INFO
db['lists']['ens'] = []
db['lists_ids']['ens'] = []
for i in db['ens']['g']:
		db['lists']['ens'].append(i)
	if i not in db['genes']:
		db['array_gene'].append(i)
		db['lists']['ens'].append(i)
		db['lists_ids']['ens'].append(len(db['array_gene'])-1)
		db['genes'][i] = {}
		db['genes'][i]['reference'] = 'ens'
		db['genes'][i]['length'] = db['ens']['g'][i]['length']
###GENELIST:BOTH
db['lists']['both'] = []
db['lists_ids']['both'] = []
for i in db['ens']['g']:
	if i in db['lists']['ref']:
		db['lists']['both'].append(i)
		db['lists_ids']['both'].append(db['array_gene'].index(i))
###GENELIST:only_ENSE
db['lists']['ens_only'] = []
db['lists_ids']['ens_only'] = []
for i in db['ens']['g']:
	if i not in db['lists']['ref']:
		db['lists']['ens_only'].append(i)
		db['lists_ids']['ens_only'].append(db['array_gene'].index(i))
###GENELIST:only_ENSE_coding_lincRNA
db['lists']['ens_only_coding'] = []
db['lists']['ens_only_lincRNA'] = []
db['lists_ids']['ens_only_coding'] = []
db['lists_ids']['ens_only_lincRNA'] = []
for i in db['ens']['g']:
	if i not in db['lists']['ref']:
		if db['ens']['g'][i]['class'] == 'protein_coding':
			db['lists']['ens_only_coding'].append(i)
			db['lists_ids']['ens_only_coding'].append(db['array_gene'].index(i))
		elif db['ens']['g'][i]['class'] == 'lincRNA':
			db['lists']['ens_only_lincRNA'].append(i)
			db['lists_ids']['ens_only_lincRNA'].append(db['array_gene'].index(i))
pickle.dump(db,open('Ref.mm9_Ref_Ens_1201.dat','wb'),True)