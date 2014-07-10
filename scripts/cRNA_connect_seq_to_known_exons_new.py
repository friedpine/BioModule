from __future__ import division
import sys
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import pairwise2
from pairwise2 import format_alignment
import subprocess
import string
import os,re,time
import MySQLdb as mb

DB_NAME = sys.argv[1]
table1 = sys.argv[2]
table2 = sys.argv[3]
species_type = sys.argv[4]
threads = sys.argv[5]
order = sys.argv[6]

def complementary(seq):
	rc = ''
	for i in seq:
		if i == 'A':
			rc=rc+'T'
		elif i == 'T':
			rc=rc+'A'
		elif i == 'C':
			rc=rc+'G'
		elif i == 'G':
			rc=rc+'C'
		else:
			rc=rc+i
	return rc

conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
cursor = conn.cursor()
#cursor.execute("select id from raw1 where self_gencode = 0 and mate_gencode = 1")
cursor.execute("select id from raw1 where self_gencode = 0")
results = cursor.fetchall()
print len(results)
ids = []
for i in results:
	#print i,int(i[0])%int(threads),order
	if int(i[0])%int(threads) == int(order):
		ids.append(int(i[0]))

for id in ids: 
	cursor.execute("select * from "+table1+" where id = %s",[id])
	results = cursor.fetchall()
	t = results[0][0:6]
	cursor.execute("select transc,ex_order from "+species_type+" where chr = %s and left_pos < %s and left_pos > %s", [t[2],t[4],t[4]-50])
	left_transcs =  cursor.fetchall()
	if left_transcs == ():
		continue 
	cursor.execute("select transc,ex_order from "+species_type+" where chr = %s and right_pos > %s and right_pos < %s", [t[2],t[5],t[5]+50])
	right_transcs = cursor.fetchall()
	if right_transcs == ():
		continue
	cursor.execute("select chr,start,start from realign where rec = 'M_G' and id = %s and chr = %s and rec = 'M_G'",[id,t[2]])
	mate_genome = tuple(set(cursor.fetchall()))
	mate_transcs = ()
	for value in mate_genome:
		cursor.execute("select transc,ex_order from "+species_type+" where chr = %s and left_pos < %s and right_pos > %s", value)
		mate_transcs +=  cursor.fetchall()
	if mate_transcs == ():
		continue
	temp = {'r':{},'m':{}}
	consensus = 'NA'
	for i in right_transcs:
		temp['r'][i[0]] = i[1]
	for i in mate_transcs:
		temp['m'][i[0]] = i[1] 
	for i in left_transcs:
		if i[0] in temp['r'] and i[0] in temp['m']:
			consensus=i[0]+'/'+str(i[1])+'/'+str(temp['r'][i[0]])+'/'+str(temp['m'][i[0]])
			cursor.execute("select * from "+species_type+" where transc = %s and ex_order = %s",[i[0],i[1]])
			t1 = cursor.fetchall()[0]
			event = t1[2]+':'+str(t1[4])+'-'+str(t1[5])
			length = t[4]-t1[4]
			cursor.execute("select * from "+species_type+" where transc = %s and ex_order = %s",[i[0],temp['r'][i[0]]])
			t1 = cursor.fetchall()[0]
			event = event+'_'+t1[2]+':'+str(t1[4])+'-'+str(t1[5])
			length = length+t1[5]-t[5]
			consen = (i[0],min(i[1],temp['r'][i[0]]),max(i[1],temp['r'][i[0]]),temp['m'][i[0]],length)
			break
	if consensus == 'NA':
		continue
	cursor.execute("insert into "+table2+" values(%s,%s,%s,%s,%s,%s,%s,%s)", (t[0],t[1],event)+consen)
	#print t,'#',event,'#',consensus,'#',left_transcs,'#',right_transcs,'#',mate_transcs
conn.commit()

