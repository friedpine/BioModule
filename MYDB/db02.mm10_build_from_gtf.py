import sys,time
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import infra01_pos2info as in1
import MySQLdb as mb
DB_NAME = 'mm10_ref2'
species = 'mm10'

gtf = "/data/Analysis/fanxiaoying/database/mm10/mm10.refseq.gtf"
#gtf = "test.gff"


sql_transc = """CREATE TABLE transc (
	id varchar(20) primary key,
	gene varchar(20),
	length INT,
	exon_count INT,
	chr varchar(20),
	strand varchar(1)
	)
"""
sql_codons = """CREATE TABLE codons (
	transc varchar(20) primary key,
	class varchar(5),
	left_pos INT,
	right_pos INT,
	unique key (transc,class)
	)
"""
sql_exons = """CREATE TABLE exons ( 
	id INT PRIMARY KEY AUTO_INCREMENT,
	chr varchar(20),
	strand varchar(1),
	left_pos INT,
	right_pos INT,
	sequence varchar(5000),
	unique key (chr,strand,left_pos,right_pos)
	)
"""
sql_introns = """CREATE TABLE introns ( 
	id INT PRIMARY KEY AUTO_INCREMENT,
	chr varchar(20),
	strand varchar(1),
	left_pos INT,
	right_pos INT,
	sequence varchar(5000),
	unique key (chr,strand,left_pos,right_pos)
	)
"""
sql_transc_exons = """CREATE TABLE t_ex (
	id varchar(20) primary key,
	ex_order INT,
#	ex_id INT,
	chr varchar(20),
	strand varchar(1),
	left_pos INT,
	right_pos INT,
	unique key (id,chr,strand,left_pos,right_pos)
#	primary key (id,ex_id),
	#foreign key (id) references transc (id),
	#foreign key (chr,strand,left_pos,right_pos) references exons (chr,strand,left_pos,right_pos)
	)
"""
values = {}
values['transc'] = []
values['start_codons'] = []
values['end_codons'] = []
values['exons'] = []
values['transc_exons'] = []

import re

if 'raw' in sys.argv:
	f = open(gtf)
	for line in f:
		t = re.split('\s+',line)
		transc = t[9][1:-2]
		if t[2] == 'CDS':
			continue
		elif t[2] == 'exon':
			 values['exons'].append((t[0],t[6],int(t[3]),int(t[4])))
			 values['transc'].append((transc,t[0],t[6]))
			 values['transc_exons'].append((transc,t[0],t[6],int(t[3]),int(t[4])))
		elif t[2] == 'start_codon':
			values['start_codons'].append((transc,int(t[3]),int(t[4])))
		elif t[2] == 'stop_codon':
			values['end_codons'].append((transc,int(t[3]),int(t[4])))
#print 'exon_counts',len(values['exons'])

conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
cursor = conn.cursor()

if not cursor.execute("show tables like 'transc'"):
	cursor.execute(sql_transc)
	cursor.executemany("""insert ignore into transc values(%s,NULL,1000,NULL,%s,%s) """,values['transc']);
	conn.commit()
if not cursor.execute("show tables like 'exons'"):
	cursor.execute(sql_exons)
	cursor.executemany("""insert ignore into exons values(NULL,%s,%s,%s,%s,NULL) """,values['exons']);
	conn.commit()
if not cursor.execute("show tables like 'codons'"):
	cursor.execute(sql_codons)
	cursor.executemany("""insert ignore into codons values(%s,'START',%s,%s) """,values['start_codons']);
	cursor.executemany("""insert ignore into codons values(%s,'STOP',%s,%s) """,values['start_codons']);
	conn.commit()
if not cursor.execute("show tables like 't_ex'"):
	cursor.execute(sql_transc_exons)
	cursor.executemany("""insert ignore into t_ex values(%s,NULL,%s,%s,%s,%s) """,values['transc_exons']);
	conn.commit()

#GET TRANSC IDS
cursor.execute("select id,strand from transc")
results = cursor.fetchall() 
transc_infos = [[i[0],i[1]] for i in results]

#BUILD_EXON_ORDER_infos
if 'BUILD_EXONS_OROER' in sys.argv:
	for transcinfo in transc_infos:
		transc = transcinfo[0]
		strand = transcinfo[1]
		if strand == '+':
			cursor.execute("select * from t_ex where id = %s order by left_pos",[transc])
		elif strand == '-':
			cursor.execute("select * from t_ex where id = %s order by right_pos DESC", [transc]) 
		results = cursor.fetchall()
		for i,j in enumerate(results):	
			cursor.execute("update t_ex set ex_order = %s where t_ex.id = %s and t_ex.chr = %s and t_ex.strand = %s and t_ex.left_pos = %s and t_ex.right_pos = %s", [i+1,j[0],j[2],j[3],j[4],j[5]])
#GET_SEQUENCE_DATA
if 'GET_SEQUENCE' in sys.argv:
	cursor.execute("select id,chr,strand,left_pos,right_pos from exons where sequence is null")
	results = cursor.fetchall()
	ranges = {}
	for i in results:
		ranges[i[0]] = i[1:]
	out = in1.get_sequenecs_from_genome_s(species,ranges)
	for i in out:
		if len(out[i])>5000:
			out[i] = out[i][0:2499]+out[i][len(out[i])-2500:len(out[i])]
		cursor.execute("update exons set sequence = %s where id = %s",[out[i],i])

conn.commit()
cursor.close()
conn.close()
