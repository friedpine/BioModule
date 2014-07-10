import sys,time
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import infra01_pos2info as in1
import MySQLdb as mb
DB_NAME = 'mm10_enseall'
species = 'mm10'

###drop tables codons,exons,g_t,gene,t_ex,transc;

gtf = "/data/Analysis/fanxiaoying/database/mm10/03.ensembl/Mus_musculus.GRCm38.73.gtf"
#gtf = "test.gff"
sql_genes = """CREATE TABLE gene (
	gene varchar(20),
	ENSG varchar(20) primary key, 
	cate varchar(30)
	)
"""	
sql_g_t = """CREATE TABLE g_t (
	ENSG varchar(20),
	transc varchar(20),
	primary key(ENSG,transc)
	) 
"""
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
	transc varchar(20),
	class varchar(5),
	left_pos INT,
	right_pos INT,
	primary key (transc,class)
	)
"""
sql_exons = """CREATE TABLE exons ( 
	id INT PRIMARY KEY AUTO_INCREMENT,
	chr varchar(20),
	strand varchar(1),
	left_pos INT,
	right_pos INT,
	length INT,
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
	transc varchar(20),
	ex_order INT,
#	ex_id INT,
	chr varchar(20),
	strand varchar(1),
	left_pos INT,
	right_pos INT,
	primary key (transc,chr,strand,left_pos,right_pos)
	)
"""
sql_transc_introns = """CREATE TABLE t_in (
	transc varchar(20),
	ex_order INT,
#	in_id INT,
	chr varchar(20),
	strand varchar(1),
	left_pos INT,
	right_pos INT,
	primary key (transc,chr,strand,left_pos,right_pos)
	)
"""

types_considered = ['protein_coding','lincRNA']

values = {}
values['genes'] = []
values['transc'] = []
values['gene_transc'] = []
values['start_codons'] = []
values['end_codons'] = []
values['exons'] = []
values['transc_exons'] = []

import re

if 'raw' in sys.argv:
	f = open(gtf)
	for line in f:
		t = re.split('\s+',line)
#		if t[2] == 'CDS' or t[1] not in types_considered:
		if t[2] == 'CDS':
			continue
		ensg = t[9][1:-2]
		transc = t[11][1:-2]
		gene = t[15][1:-2]
		t[0] = 'chr'+t[0]
		if t[2] == 'exon':
			 values['genes'].append((gene,ensg,t[1]))
			 values['gene_transc'].append((ensg,transc))
			 values['exons'].append((t[0],t[6],int(t[3]),int(t[4]),int(t[4])-int(t[3])+1))
			 values['transc'].append((transc,t[0],t[6]))
			 values['transc_exons'].append((transc,int(t[13][1:-2]),t[0],t[6],int(t[3]),int(t[4])))
		if t[2] == 'start_codon':
			values['start_codons'].append((transc,int(t[3]),int(t[4])))
		if t[2] == 'stop_codon':
			values['end_codons'].append((transc,int(t[3]),int(t[4])))
#print 'exon_counts',len(values['exons'])

conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
cursor = conn.cursor()

if not cursor.execute("show tables like 'gene'"):
	cursor.execute(sql_genes)
	cursor.executemany("""insert ignore into gene values(%s,%s,%s) """,values['genes']);
	conn.commit()
if not cursor.execute("show tables like 'g_t'"):
	cursor.execute(sql_g_t)
	cursor.executemany("""insert ignore into g_t values(%s,%s) """,values['gene_transc']);
	conn.commit()
if not cursor.execute("show tables like 'transc'"):
	cursor.execute(sql_transc)
	cursor.executemany("""insert ignore into transc values(%s,NULL,1000,NULL,%s,%s) """,values['transc']);
	conn.commit()
if not cursor.execute("show tables like 'exons'"):
	cursor.execute(sql_exons)
	cursor.executemany("""insert ignore into exons values(NULL,%s,%s,%s,%s,%s,NULL) """,values['exons']);
	conn.commit()
if not cursor.execute("show tables like 'codons'"):
	cursor.execute(sql_codons)
	cursor.executemany("""insert ignore into codons values(%s,'START',%s,%s) """,values['start_codons']);
	cursor.executemany("""insert ignore into codons values(%s,'STOP',%s,%s) """,values['start_codons']);
	conn.commit()
if not cursor.execute("show tables like 't_ex'"):
	cursor.execute(sql_transc_exons)
	cursor.executemany("""insert ignore into t_ex values(%s,%s,%s,%s,%s,%s) """,values['transc_exons']);
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
			cursor.execute("select * from t_ex where transc = %s order by left_pos",[transc])
		elif strand == '-':
			cursor.execute("select * from t_ex where transc = %s order by right_pos DESC", [transc]) 
		results = cursor.fetchall()
		for i,j in enumerate(results):	
			cursor.execute("update t_ex set ex_order = %s where t_ex.transc = %s and t_ex.chr = %s and t_ex.strand = %s and t_ex.left_pos = %s and t_ex.right_pos = %s", [i+1,j[0],j[2],j[3],j[4],j[5]])
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
if 'BUILD_INTRONS' in sys.argv:
	values = []
	cursor.execute("SELECT * FROM t_ex ORDER BY transc,ex_order")
	r0 =  cursor.fetchall()
	print len(r0)
	for i in range(0,len(r0)-1):
		if r0[i][0] == r0[i+1][0] and r0[i][1]+1 == r0[i+1][1]:
			if r0[i][3] == '-':
				values.append([r0[i][0],r0[i][1],r0[i][2],r0[i][3],r0[i+1][5]+1,r0[i][4]-1])
			if r0[i][3] == '+':
				values.append([r0[i][0],r0[i][1],r0[i][2],r0[i][3],r0[i][5]+1,r0[i+1][4]-1])
	print len(values)
	if not cursor.execute("show tables like 't_in'"):
		cursor.execute(sql_transc_introns)
	cursor.executemany("insert into t_in values(%s,%s,%s,%s,%s,%s)",values)
if 'exon_counts' in sys.argv:
	cursor.execute("SELECT COUNT(*),sum(left_pos),sum(right_pos),transc FROM t_ex a GROUP BY transc")
	r0 = cursor.fetchall()
	values = []
	for i in r0:
		values.append([i[0],int(i[2])-int(i[1])+int(i[0]),i[3]])
	print values[1:10]
	cursor.executemany("update transc set exon_count = %s,length=%s where id = %s",values)

conn.commit()
cursor.close()
conn.close()
