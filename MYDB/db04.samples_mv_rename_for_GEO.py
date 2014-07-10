import MySQLdb as mb
import subprocess,re,sys

conn=mb.connect(host="localhost",user="root",passwd="123456",db='p00_chips')
cursor = conn.cursor()

if 'rename_samples' in sys.argv:
	f = open('/datb/fanxiaoying/project/project04_RNA_chip/GEO/samples')
	values = []
	for line in f:
		a = re.split('\s+',line)
		values.append([a[2],a[0]])
	cursor.executemany("update samples set publish2 = %s where sample = %s",values)
	conn.commit()

cursor.execute("SELECT * FROM p00_chips.samples WHERE publish2 != ''")
r0 = cursor.fetchall()
pathes = []	

if 'mv_fastq_files' in sys.argv:
	for i in r0:
		cursor.execute("SELECT path from files where sample = %s and type=%s",[i[0],'fq1'])
		path0 = cursor.fetchall()[0][0]
		path1 = "/datb/fanxiaoying/project/project04_RNA_chip/GEO/friedpine@gmail.com/sample_"+i[3]+"_R1.fastq.gz"
		#subprocess.call("mv "+path0+" "+path1,shell='True')
		pathes.append([path1,path0])
		print 'mv',path0,path1
		cursor.execute("SELECT path from files where sample = %s and type=%s",[i[0],'fq2'])
		path0 = cursor.fetchall()[0][0]
		path1 = "/datb/fanxiaoying/project/project04_RNA_chip/GEO/friedpine@gmail.com/sample_"+i[3]+"_R2.fastq.gz"
		#subprocess.call("mv "+path0+" "+path1,shell='True')
		pathes.append([path1,path0])
		print 'mv',path0,path1
	cursor.executemany("update files set path = %s where path = %s",pathes)
	conn.commit()
if 'cp_fpkm_files_for_shen' in sys.argv:
	f = open("/datb/fanxiaoying/project/project04_RNA_chip/GEO/good_cells")
	for line in f:
		r0 = re.split('\s+',line)[0]
		cursor.execute("SELECT path from files where type = 'cuff_ensembl' and sample = %s",[r0])
		path0 = cursor.fetchall()[0][0]
		cursor.execute("SELECT publish2 from samples where sample = %s",[r0])
		path1= "/datb/fanxiaoying/project/project04_RNA_chip/GEO/for_shen/FPKM_"+cursor.fetchall()[0][0]+'.txt'
		print "cp "+path0+' '+path1

if 'info_genes_rpkm' in sys.argv:
	for i in r0:
		print i[0],i[2]