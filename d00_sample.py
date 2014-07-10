import MySQLdb as mb
def get_path(dbname,sample,type):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=dbname)
	cursor = conn.cursor()
	cursor.execute("select path from files where sample = %s and type = %s",([sample,type]))
	return cursor.fetchall()[0][0]
def cp_move_files(pathes,operation,rec):
	if operation == 'cp':
		for i in pathes:
			subprocess.call("cp "+i[0]+" "+i[1],shell='True')
	if operation == 'mv':
		for i in pathes:
			subprocess.call("mv "+i[0]+" "+i[1],shell='True')
def cufflinks_command(dbname,samples,bamrec,outrec,outdirrec,gtf):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=dbname)
	cursor = conn.cursor()
	out = []
	for sample in samples:
		cursor.execute("select path from files where sample = %s and type = %s",([sample,bamrec]))
		bam = cursor.fetchall()[0][0]
		cursor.execute("select path from files where sample = %s and type = %s",([sample,outdirrec]))
		outdir = cursor.fetchall()[0][0]
		cursor.execute("insert ignore into files values(%s,%s,%s,%s)",([sample,outrec,outdir+'/genes.fpkm_tracking',8888]))
		conn.commit()
		out.append('/data/Analysis/fanxiaoying/software/cufflinks-2.1.1.Linux_x86_64/cufflinks -p 4 -o '+outdir+' -G '+gtf+' '+bam)
	return out
def pair_end_insert_size(dbname,samples,bamrec,outdir,outrec):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=dbname)
	cursor = conn.cursor()
	out = []
	for sample in samples:
		cursor.execute("select path from files where sample = %s and type = %s",([sample,bamrec]))
		bam = cursor.fetchall()[0][0]
		outfile = outdir+'/insert_size_'+sample+'.txt'
		cursor.execute("insert ignore into files values(%s,%s,%s,%s)",([sample,outrec,outfile,8888]))
		conn.commit()
		out.append('samtools view -f 2 '+bam+" | awk '{print $9}' > "+outdir+'/insert_size_'+sample+'.txt')
	return out
def get_sample_file(cursor,sample,type):
	cursor.execute("select path from files where sample = %s and type = %s",[sample,type])
	return cursor.fetchall()[0][0]
def get_sample_info(cursor,sample,type):
	cursor.execute("select "+type+" from samples where sample = %s",[sample])
	return cursor.fetchall()[0][0]
def insert_sample_file(cursor,sample,type,path):
	cursor.execute("insert into files values(%s,%s,%s,NULL)",[sample,type,path])
def table_2_dict(cursor,tablename,columes):
	gt = {}
	cursor.execute("select %s,%s from %s" %(columes[0],columes[1],tablename))
	r0 = cursor.fetchall()
	for i in r0:
		gt[i[0]] = i[1]
	return gt
def run_pipeline(cmdname,n):
	import subprocess,time
	handle = []
	for i in range(0,n):
		handle.append('PP'+str(i))
	info = cmdname
	if len(info)>n:
		next = n
		for i in range(0,n):
			 handle[i]=subprocess.Popen(info[i],shell='True')
		while 1:
			running_count = 0
			for i in range(0,n):
				if handle[i].poll() is None:
					print str(i)+" Is RUNNING"
					running_count += 1
				elif next < len(info):
					handle[i] = subprocess.Popen(info[next],shell='True')
					print str(i)+"NEW_tasks for this handle"
					next += 1
				time.sleep(1)
			if running_count == 0:
				break
	else:
		for i in range(0,len(info)):
			handle[i]=subprocess.Popen(info[i],shell='True')
		while 1:
			running_count = 0
			for i in range(0,len(info)):
				if handle[i].poll() is None:
					print str(i)+" Is RUNNING"
					running_count += 1
					time.sleep(1)
			if running_count == 0:
				br