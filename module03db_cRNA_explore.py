import sys,os
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import random as rd
import numpy as np
#import module03_rna_procerss as m03
import infra00_ranges_operate as in0
import infra01_pos2info as in1
import infra02_blast_info_proc as in2
import infra03_conservation as in3
import MySQLdb as mb
import d02_db_tables as d02

def random_build_fake_cRNA(refdb,dbname,table1,numbers):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=dbname)
	cursor = conn.cursor()
	cursor.execute("SELECT id,exon_count FROM "+refdb+".transc WHERE exon_count>5 AND id NOT IN (SELECT transc FROM mm_cRNA.01_all_events)")
	r0 = list(cursor.fetchall())
	transc_3k = rd.sample(r0,numbers)
	values = []
	for transc in transc_3k:
		r1 = rd.sample(range(2,transc[1]),2)
		values.append([transc[0],min(r1),max(r1)])
	cursor.executemany("insert into "+table1+" (id,transc,ex_left_id,ex_right_id)values(NULL,%s,%s,%s)",values)
	conn.commit()

def blast_introns_sequence(refdb,table_in,dbname,table1,table2,evalue,wordsize):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=dbname)
	cursor = conn.cursor()
	cursor.execute("select * from "+table1)
	r0 = cursor.fetchall()
	for t in r0:
		print '###########',t
		try:
			ranges = {}
			cursor.execute("select * from "+refdb+'.'+table_in+" where transc = %s and ex_order = %s ",[t[1],t[2]-1])
			r1 = cursor.fetchall()[0]
			ranges['r2'] = [r1[2],'+',r1[4],r1[5]]
			cursor.execute("select * from "+refdb+'.'+table_in+" where transc = %s and ex_order = %s ",[t[1],t[3]])
			r1 = cursor.fetchall()[0]
			ranges['r1'] = [r1[2],'+',r1[4],r1[5]]
			in2.blast_genome_positions('mm10',ranges,evalue,wordsize,'RC')
			result = in2.blast_fmt7_out_read_for_db('./temp_blast_r1r2.txt',40,'RC')
			if result != []:
				for r2 in result:
					cursor.execute("insert into "+table2+" values(%s,%s,%s,%s,%s,%s,%s,%s,%s)",[t[0]]+r2)
		except:
			cursor.execute("insert into "+table2+" values(%s,0,0,0,0,0,0,0,0)",[t[0]])
		conn.commit()

def get_last_junction_reads_counts(cursor,conn,table_events,table_last_junc,table_files,file_format):
	check_col_1 = d02.check_table_columes(cursor,table_events,['id','sample','transc'])
	if check_col_1 != []:
		print "THESE_COLUMES_DOSENT_EXISTS"," ".join(check_col_1)
	check_col_2 = d02.check_table_columes(cursor,table_events,['transc','chr','left_pos','right_pos'])
	if check_col_2 != []:
		print "THESE_COLUMES_DOSENT_EXISTS"," ".join(check_col_2)
	sql = "select  a.id,a.sample,b.chr,b.left_pos,b.right_pos from "+table_events+"a join "+table_last_junc+" b on a.transc = b.transc"
	cursor.execute(sql)
	events = cursor.fetchall()
	events_dict = {}
	for x in events:
		if x[1] not in events_dict:
			events_dict[x[1]] = []
		events_dict[x[1]].append()
	samples = events_dict.keys()
	for sample in samples:
