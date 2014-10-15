from __future__ import division
import module02_coverage as m02
import subprocess
import cPickle as pickle
import time
import string
import matplotlib.pyplot as plt

###FUNC:build_bulk_coverage(self,counts_file,bamfile)
###FUNC:normalize_by_bulk(self,bulk,sample):
###FUNC:bin_by_length(self,sample,geneslist,normal_way)
###FUNC:plot_all_six_lines(self,samples,labels,legends,colors,record)
try:
	t = pickle.load(open("Data.coverage_normalize.dat","rb"))
except:
	db = {}
	db = pickle.load(open("mm9_Ref_Ens.dat","rb"))
	print db['ref'].keys()
	t = m02.coverage()
	t.initialize(db)
#	t.build_coverage('bulk','sum.G01.txt','00.data/G01.sort.bam')
	t.build_coverage('S01','sum.B13.txt','00.data/B13.sort.bam')
	t.build_coverage('S54','Sum.S54.txt','00.data/S54.bam')
	t.bin_by_length('S01','well_expressed_list','raw')
	t.bin_by_length('S54','well_expressed_list','raw')
	pickle.dump(t,open("Data.coverage_normalize.dat","wb"),True)

labels = ['0-1K','1-2K','2-3K','3-4K','4-5K','>5K']
legends = ['NEW','TANG']
colors = ['b','r']

t.plot_all_six_lines(['S01',"S54"],labels,legends,colors,'depth_10_single_isoform')

f = open('info2.txt','w')
print >>f,"Bin by length"
for i in t['S01']['bin']:
	print >>f,'mean',i
for i in t['S01']['bin_std']:
        print >>f,"std",i
f.close


#20140824 Today is my birthday! I am now 24 years old! How time flies!
