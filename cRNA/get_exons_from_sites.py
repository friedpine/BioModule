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

sites_db = sys.argv[1]
sites_table = sys.argv[2]
species_type = sys.argv[3]
species_type_in = sys.argv[4]

chr_p = int(sys.argv[5])
pos_p = int(sys.argv[6])

conn=mb.connect(host="localhost",user="root",passwd="123456",db=sites_db)
cursor = conn.cursor()
try:
	cursor.execute("ALTER TABLE "+sites_table+" ADD transc VARCHAR(20),ADD ex_l INT,ADD ex_r INT,ADD annotype VARCHAR(10)")
except:
	print "EXISTS"

cursor.execute("select * from "+sites_table)
results = cursor.fetchall()
for r0 in results:
	t = r0
	print t,
	cursor.execute("select transc,ex_order from "+species_type+" where chr = %s and left_pos >=%s and left_pos <=%s", [t[chr_p],t[chr_p+1]-2,t[chr_p+1]+2])
	left_transcs =  cursor.fetchall()
	left_type = "S"
	if left_transcs == ():
		cursor.execute("select transc,ex_order from "+species_type+" where chr = %s and left_pos <%s and right_pos>%s", [t[chr_p],t[chr_p+1],t[chr_p+1]])
		left_transcs =  cursor.fetchall()
		left_type = "Ex"
		if left_transcs == ():
			cursor.execute("select transc,ex_order from "+species_type_in+" where chr = %s and left_pos <=%s and right_pos>=%s", [t[chr_p],t[chr_p+1],t[chr_p+1]])
			left_transcs =  cursor.fetchall()
			left_type = "In"
			if left_transcs == ():
				left_type = "N"
	print left_transcs, 
	cursor.execute("select transc,ex_order from "+species_type+" where chr = %s and right_pos >=%s and right_pos <=%s", [t[chr_p],t[chr_p+2]-2,t[chr_p+2]+2])
	right_transcs = cursor.fetchall()
	right_type = "S"
	if right_transcs == ():
		cursor.execute("select transc,ex_order from "+species_type+" where chr = %s and left_pos <%s and right_pos>%s", [t[chr_p],t[chr_p+2],t[chr_p+2]])
		right_transcs =  cursor.fetchall()
		right_type = "Ex"
		if right_transcs == ():
			cursor.execute("select transc,ex_order from "+species_type_in+" where chr = %s and left_pos <=%s and right_pos>=%s", [t[chr_p],t[chr_p+2],t[chr_p+2]])
			right_transcs =  cursor.fetchall()
			right_type = "In"
			if right_transcs == ():
				right_type = "N"
	print right_transcs,
	temp = {'r':{}}
	consensus = 'NA'
	for i in right_transcs:
		temp['r'][i[0]] = i[1]
	for i in left_transcs:
		if i[0] in temp['r']:
			cursor.execute("update "+sites_table+" set transc=%s,ex_l=%s,ex_r=%s where pos= %s",[i[0],str(i[1]),str(temp['r'][i[0]]),t[pos_p]])
			print "OK NICE!!"
			break
	cursor.execute("update "+sites_table+" set annotype=%s where pos= %s",[left_type+"_"+right_type,t[pos_p]])
	conn.commit()



