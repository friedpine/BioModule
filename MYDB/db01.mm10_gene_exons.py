import sys
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import infra01_pos2info as in1
import cPickle as pickle
db = pickle.load(open('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/Ref.exon.mm10.dat'))
import MySQLdb as mb
DB_NAME = 'TEST01'


conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
cursor = conn.cursor()

sql_transc = """CREATE TABLE transc (
	id varchar(20) PRIMARY KEY,
	gene varchar(20),
	length INT,
	exon_count INT,
	chr varchar(10),
	starnd varchar(1),
	left_pos INT,
	right_pos INT
	)
"""
sql_exons = """CREATE TABLE exons (
	id INT PRIMARY KEY AUTO_INCREMENT,
	chr varchar(10),
	starnd varchar(1),
	left_pos INT,
	right_pos INT,
	sequence varchar(5000))
"""
if not cursor.execute("show tables like 'transc'"):
	values = []
	for t in db['transc_info']:
		info = db['transc_info'][t]
		if len(info['exons_pos'])>0:
			values.append((t,'UNSET',1000,info['count'],info['chr'],info['strand'],info['exons_pos'][0],info['exons_pos'][-1]))
		else:
			print info
	cursor.execute(sql_transc)
	cursor.executemany("""insert into transc values(%s,%s,%s,%s,%s,%s,%s,%s) """,values);
if not cursor.execute("show tables like 'exons'"):
	values = []
	for chr in db['exon_info']:
		for left in db['exon_info'][chr]:
			for right in db['exon_info'][chr][left]:
				info = db['exon_info'][chr][left][right]
				if 'strand' in info and 'seq' in info:
					values.append((chr,info['strand'],left,right,info['seq']))
				else:
					print info
	cursor.execute(sql_exons)
	cursor.executemany("""insert into exons values(NULL,%s,%s,%s,%s,%s) """,values);

conn.commit()
cursor.close()
conn.close()



