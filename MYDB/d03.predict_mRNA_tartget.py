import sys,time,re
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import infra01_pos2info as in1
import MySQLdb as mb
DB_NAME = 'mm10_ensembl'

def reverse_comp(seq):
	out = ''
	for i in seq[::-1]:
		if i == 'A':
			out += 'T'
		elif i == 'U' or i=='T':
			out += 'A'
		elif i == 'C':
			out += 'G'
		elif i == 'G':
			out += 'C'
	return out
#print reverse_comp('UCACAUU')

conn=mb.connect(host="localhost",user="root",passwd="123456",db='mm10_ensembl')
cursor = conn.cursor()

if sys.argv[1] == 'predict_exon_mRNA':
	f = open('/tmp/7mers.txt')
	for line in f:
		miRNAid = int(re.split('\s+',line)[0])
		if miRNAid%int(sys.argv[1]) != int(sys.argv[2]):
			continue
		seq7 = reverse_comp(re.split('\s+',line)[1])
		cursor.execute("select id,sequence from exons where locate(%s,sequence)>0 ",[seq7])
		r = cursor.fetchall()
		values = []
		for i in r:
			t = re.findall(seq7,i[1])
			values.append((i[0],miRNAid,len(t)))
		cursor.executemany("insert into ex_miRNA_count values(%s,%s,%s)",values)
		conn.commit()
if sys.argv[1] == 'count_ex_sites':
	cursor.execute("select ex,sum(count) from ex_miRNA_count group by ex")
	r = cursor.fetchall()
	values = []
	for i in r:
		values.append((i[1],i[0]))
	cursor.executemany("update exons set miRNA_count = %s where id = %s ",values)
	conn.commit()
if sys.argv[1] == 'transc_miRNA_estimate':
	db_ref = 'mm10_ensembl'
	values = []
	cursor.execute("select id from transc")
	r1 = cursor.fetchall()
	for index,i in enumerate(r1):
		if index%int(sys.argv[2]) == int(sys.argv[3]):
			cursor.execute("select * from "+db_ref+".t_ex where transc = %s",[i[0]])
			r2 = cursor.fetchall()
			length = 0
			sites = 0
			for j in r2:
				length += j[5]-j[4]
				cursor.execute("select miRNA_count from "+db_ref+".exons where chr=%s and left_pos = %s and right_pos = %s",[j[2],j[4],j[5]])
				sites += int(cursor.fetchall()[0][0])
			value = [int(length),int(sites),i[0]]
			values.append((int(length),sites,i[0]))
	#print values
	cursor.executemany("update transc set length = %s ,miRNA_sites = %s where id = %s ",values)
	conn.commit()
	# cursor.close()
	# conn.close()
