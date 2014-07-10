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
anchor_len = sys.argv[5]
threads = sys.argv[6]
order = sys.argv[7]

chr_prefix = ''
if re.match("mm",species_type):
	chr_prefix = ''

sql = "CREATE TABLE " +table2+"""(
  `id` int(11) DEFAULT NULL,
  `sample` varchar(100) DEFAULT NULL,
  `event` varchar(100) DEFAULT NULL,
  `transc` varchar(100) DEFAULT NULL,
  `l_ex` int(11) DEFAULT NULL,
  `r_ex` int(11) DEFAULT NULL,
  `m_ex` int(11) DEFAULT NULL,
  `l_dis` int(11) DEFAULT NULL,
  `r_dis` int(11) DEFAULT NULL,
  `t_dis` int(11) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1"""

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

try:
	cursor.execute(sql)
except:
	print "EXISTS"

cursor.execute("select id from "+table1+" where self_gencode = 0")
results = cursor.fetchall()
print len(results)
ids = []
for i in results:
	if int(i[0])%int(threads) == int(order):
		ids.append(int(i[0]))
for id in ids: 
	cursor.execute("select * from "+table1+" where id = %s",[id])
	results = cursor.fetchall()
	t = list(results[0][0:7])
	t[2] = chr_prefix+str(t[2])
	cursor.execute("select transc,ex_order from "+species_type+" where chr = %s and left_pos < %s and right_pos > %s", [t[2],t[4]+2,t[4]-2])
	left_transcs =  cursor.fetchall()
	print left_transcs
	if left_transcs == ():
		continue 
	cursor.execute("select transc,ex_order from "+species_type+" where chr = %s and left_pos < %s and right_pos > %s", [t[2],t[5]+2,t[5]-2])
	right_transcs = cursor.fetchall()
	if right_transcs == ():
		continue
	temp = {'r':{}}
	consensus = 'NA'
	for i in right_transcs:
		temp['r'][i[0]] = i[1]
	info = []
	for i in left_transcs:
		if i[0] in temp['r']:
			consensus=i[0]+'/'+str(i[1])+'/'+str(temp['r'][i[0]])
			cursor.execute("select * from "+species_type+" where transc = %s and ex_order = %s",[i[0],i[1]])
			t1 = cursor.fetchall()[0]
			event = t1[2]+':'+str(t1[4])+'-'+str(t1[5])
			length1 = t[4]-t1[4]
			cursor.execute("select * from "+species_type+" where transc = %s and ex_order = %s",[i[0],temp['r'][i[0]]])
			t1 = cursor.fetchall()[0]
			length2 = t1[5]-t[5]
			if info == []:
				info = [i[0],min(i[1],temp['r'][i[0]]),max(i[1],temp['r'][i[0]]),0,length1,length2,length1+length2+2*int(anchor_len)-len(t[6])]
				event_s = event+'_'+t1[2]+':'+str(t1[4])+'-'+str(t1[5])
			elif abs(length1+length2+2*int(anchor_len)-len(t[6]))<abs(info[-1]):
				info = [i[0],min(i[1],temp['r'][i[0]]),max(i[1],temp['r'][i[0]]),0,length1,length2,length1+length2+2*int(anchor_len)-len(t[6])]
				event_s = event+'_'+t1[2]+':'+str(t1[4])+'-'+str(t1[5])
	if consensus == 'NA':
		continue
	cursor.execute("insert into "+table2+" values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)", [t[0],t[1],event_s]+info)
	conn.commit()


