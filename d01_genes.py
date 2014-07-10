import os,re
def fpkm_of_genes(refdb,reftable,refcolum,project_db,samples,cuffrec,table1):
	import MySQLdb as mb
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=project_db)
	cursor = conn.cursor()
	cursor.execute("show tables like %s",[table1])
	r0 = cursor.fetchall()
	if r0 == ():
		cursor.execute("create table "+table1+" (id varchar(30) primary key)")
		conn.commit()
		cursor.execute("select "+refcolum+" from "+refdb+'.'+reftable)
		r0 = cursor.fetchall()
		transcs = [[i[0]] for i in r0]
		cursor.executemany("insert ignore into "+table1+" values(%s)",transcs)
		conn.commit()
	for sample in samples:
		cursor.execute("explain "+project_db+"."+table1)
		r0 = cursor.fetchall()
		table_columes = [i[0] for i in r0]
		print table_columes
		sample = 'mm_'+re.sub('-','_',sample)
		if sample not in table_columes:
			print sample
			cursor.execute("alter table "+project_db+'.'+table1+" add "+str(sample)+" float default 0")
			conn.commit()
	for sample in samples:
		values = []
		cursor.execute("select path from files where type = %s and sample = %s",[cuffrec,sample])
		cuff_file = cursor.fetchall()[0][0] 
		if not os.path.exists(cuff_file):
			print sample,"CUFF_FILE_NOT_EXISTS"
			continue
		f = open(cuff_file)
		f.readline()
		for line in f:
			t = re.split('\s+',line)
			if float(t[9]) > 0:
				values.append([float(t[9]),t[0]])
		sample = 'mm_'+re.sub('-','_',sample)
		cursor.executemany("update "+table1+" set "+sample+"= %s where id = %s",values)
		conn.commit()