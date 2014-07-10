import re, sys, os, copy
import subprocess
import time
import cPickle as pickle
import MySQLdb as mb
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import infra00_ranges_operate as in0


DB_NAME = sys.argv[1]
table_exons = sys.argv[2]
table_t_ex = sys.argv[3]
table_transc = sys.argv[4]
table_ref = sys.argv[5]


conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
cursor = conn.cursor()
exons = []
t_ex = []
transc = []

if "STEP1" in sys.argv:
	psl = sys.argv[6]
	file = open(psl)
	for i in range(5):
		file.readline()
	for line in file:
		t = re.split('\s+',line)
		if len(t)<10 or int(t[0])<0.95*int(t[10]):
			continue
		chr = t[13]
		exon_count = int(t[17])
		exons_start = [int(i) for i in re.split(',',t[20])[:-1]]
		exons_size = [int(i) for i in re.split(',',t[18])[:-1]]
		exons_ranges = [[exons_start[i],exons_start[i]+exons_size[i]-1] for i in range(exon_count)]
		exons_ranges = in0.merge_ranges(exons_ranges,25)
		info = ''
		exonsize = ''
		intronsize = ''
		for i in exons_ranges:
			exons.append([chr,i[0],i[1],i[1]-i[0]])
			t_ex.append([t[9],chr,i[0],i[1]])
			info += str(i[0])+'-'+str(i[1])+','
			exonsize += str(i[1]-i[0])+','
		for i in range(1,len(exons_ranges)):
			intronsize += str(exons_ranges[i][0]-exons_ranges[i-1][1]-1)+','
		transc.append([t[9],t[10],t[13],t[15],t[16],len(exons_ranges),info,exonsize,intronsize])
	cursor.executemany("insert ignore into "+table_exons+" (chr,left_pos,right_pos,length)values(%s,%s,%s,%s)",exons)
	cursor.executemany("insert ignore into "+table_t_ex+" (transc,chr,left_pos,right_pos)values(%s,%s,%s,%s)",t_ex)
	cursor.executemany("insert ignore into "+table_transc+" (id,length,chr,left_pos,right_pos,exon_count,info,exon_size,intron_size)values(%s,%s,%s,%s,%s,%s,%s,%s,%s)",transc)
	conn.commit()

##GET THE RELATIONSHIP WITH OTHER EXONS

cursor.execute("select id,chr,left_pos,right_pos from "+table_exons)
r0 = cursor.fetchall()
count = 0
infos = []
for i in r0:
	info = [0,0,0]
	count += 1
	if count%1000 == 1 or count == len(r0):
		cursor.executemany("update "+table_exons+ " set left_exon_dis = %s,right_exon_dis = %s where id = %s",infos)
		conn.commit()
		infos = []
		print count
	cursor.execute("select id,left_pos,right_pos from "+table_ref+" where chr = %s and left_pos>%s and left_pos<%s and right_pos<%s",[i[1],i[2]-100000,i[2],i[3]+5000])
	r1 = cursor.fetchall()
	if r1 == ():
		info[0] = 100000
	cursor.execute("select id,left_pos,right_pos from "+table_ref+" where chr = %s and right_pos>%s and right_pos<%s and left_pos>%s",[i[1],i[3],i[3]+100000,i[3]-5000])
	r1 = cursor.fetchall()
	if r1 == ():
		info[1] = 100000
	info[2] = i[0]
	infos.append(info)




	
