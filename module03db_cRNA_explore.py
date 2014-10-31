import sys,os
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import random as rd
import numpy as np
import infra00_ranges_operate as in0
import infra01_pos2info as in1
import infra02_blast_info_proc as in2
import infra03_conservation as in3
import infra06_bam_handle as in6
import MySQLdb as mb
import d00_sample as d00
import d02_db_tables as d02
import pysam

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

def blast_introns_sequence(cursor,conn,tablein,table_result,species,ref,evalue,wordsize,cut,folder,server="TANG"):
	try:
		cursor.execute('create table '+table_result+" as select * from mm_crna_explore.blast_template")
	except:
		print "EXISTS!"
	cursor.execute("select event,chr,in1_s,in1_e,in2_s,in2_e from "+tablein)
	r0 = cursor.fetchall()
	for t in r0:
		try:
			pos_list1 = [[t[1],t[2],t[3],'+']]
			pos_list2 = [[t[1],t[4],t[5],'+']]
			list1_name = ['up']
			list2_name = ['down']
			record = t[0]
			values = in2.blast_ref_positions(cursor,species,ref,pos_list1,pos_list2,list1_name,list2_name,evalue,wordsize,cut,folder,record,'RC',server)
		except:
			print "Failed To Get BLAST Result!"
			continue
		print values
		cursor.executemany("insert into "+table_result+" values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",values)
		conn.commit()

def get_last_junction_reads_counts(cursor,conn,table_events,table_last_junc,table_files,filetype,chr_prefix):
	check_col_1 = d02.check_table_columes(cursor,table_events,['id','sample','transc'])
	if check_col_1 != []:
		print "THESE_COLUMES_DOSENT_EXISTS"," ".join(check_col_1)
	check_col_2 = d02.check_table_columes(cursor,table_last_junc,['transc','chr','left_pos','right_pos'])
	if check_col_2 != []:
		print "THESE_COLUMES_DOSENT_EXISTS"," ".join(check_col_2)
	sql = "select  a.id,a.sample,b.chr,b.left_pos,b.right_pos from "+table_events+" a join "+table_last_junc+" b on a.transc = b.transc"
	cursor.execute(sql)
	events = cursor.fetchall()
	events_dict = {}
	for x in events:
		if x[1] not in events_dict:
			events_dict[x[1]] = []
		events_dict[x[1]].append(x)
	samples = events_dict.keys()
	if not d00.check_validness_of_bamfiles(cursor,table_files,samples,filetype):
		return 0
	time = 0
	for sample in samples:
		bampath = d00.get_path1(cursor,table_files,sample,filetype)
		bamfile = pysam.Samfile(bampath,"rb")
		ids = []
		data = []
		for event in events_dict[sample]:
			counts = in6.get_junction_reads_counts(bamfile,event[2][chr_prefix:],event[3],event[4])
			data.append(counts)
			ids.append(event[0])	
		d02.append_colume_info_to_tables(cursor,conn,table_events,'last_exon_junction','INT',ids,data)
		print "FINISHED",sample

def get_sequence_of_splicing_site(cursor,conn,species,intable,outtable):
	cursor.execute("select event,strand,chr,in1_s,in1_e,in2_s,in2_e from "+intable)
	infos = cursor.fetchall()
	d02.Create_Mysql_tables(cursor,conn,outtable,['event','up5','up3','down5','down3'],['varchar(60)','varchar(10)','varchar(10)','varchar(10)','varchar(10)'])
	out = []
	for x in infos:
		line = [x[0]]
		chr = x[2]
		strand = x[1]
		if x[1] == "-":
			poses = [[chr,x[4]-3,x[4]+4,strand],[chr,x[3]-4,x[3]+3,strand],[chr,x[6]-3,x[6]+4,strand],[chr,x[5]-4,x[5]+3,strand]]
			line=line+in1.get_pos_seqs(cursor,species,'genome',poses,"TANG")
		if x[1] == "+":
			poses = [[chr,x[3]-4,x[3]+3,strand],[chr,x[4]-3,x[4]+4,strand],[chr,x[5]-4,x[5]+3,strand],[chr,x[6]-3,x[6]+4,strand]]
			line=line+in1.get_pos_seqs(cursor,species,'genome',poses,"TANG")
		out.append(line)
	d02.Insert_lines_Into_table(cursor,conn,outtable,['event','up5','up3','down5','down3'],out)


