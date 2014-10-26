import re, sys, os
import subprocess
import time
import d00_sample as d00
import d02_db_tables as d02

def add_files(cmd,file):
	f = open(file)
	out = "#CMD\n"+cmd+"\n"+"#FILE\t"+file+"\n"
	for i in f:
		out+=i
	return out

def CUFFLINKS(cursor,conn,species,ref,samples,intype,folder,rec,server="TANG"):
	cmds = []
	for sample in samples:
		insert = rec+'_'+sample
		outdir = folder+'/'+insert
		if not os.path.exists(outdir):
			os.mkdir(outdir)
		bam = d00.get_sample_file(cursor,sample,intype)
		path = outdir+'/genes.fpkm_tracking'
		refpath = d00.get_ref(cursor,species,'gtf',ref,server)
		cmd = '/data/Analysis/fanxiaoying/software/cufflinks_221/cufflinks -o %s -p 6 -G %s %s' %(outdir,refpath,bam)
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,cmd])
		cmds.append(cmd)
	conn.commit()
	return cmds

def CUFFDIFF(cursor,conn,species,ref,G1,G2,G1name,G2name,intype,folder,rec):
	path = folder+"/"+rec
	cmd = "/data/Analysis/fanxiaoying/software/cufflinks_221/cuffdiff -p 8 -o "+path+" "+refall[species]['gtf'][ref]+" -L "+G1name+","+G2name+" "
	for sample in G1:
		cmd = cmd+d00.get_sample_file(cursor,sample,intype)+","
	cmd = cmd[:-1]+" "
	for sample in G2:
		cmd = cmd+d00.get_sample_file(cursor,sample,intype)+","
	print cmd[:-1]
	cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",["CUFFDIFF",rec,path,cmd])
	conn.commit()

def FPKM_DB(cursor,conn,genes_table,samples,intype,insert,tablename):
	sql = 'create table %s select * from %s' %(tablename,genes_table)
	genes = d00.table_2_dict(cursor,genes_table,['gene','gene'])
	try:
		cursor.execute(sql)
		conn.commit()
		cursor.execute("create index reci on "+tablename+"(gene);")
		conn.commit()
	except:
		print "exists"
	for sample in samples:
		exp = d00.get_sample_info(cursor,sample,'exp')+insert
		try:
			cursor.execute("alter table "+tablename+" add "+exp+" float DEFAULT '0'")
			conn.commit()
		except:
			print "EXISTS",
		fpkm = []
		f = open(d00.get_sample_file(cursor,sample,intype))
		f.readline()
		print d00.get_sample_file(cursor,sample,intype)
		for line in f:
			t = re.split('\s+',line)
			if re.findall('\.',t[0]):
				t[0] = re.split('\.',t[0])[0]
			if t[9] > 0 and t[0] in genes:
				fpkm.append([t[9],t[0]])
		print fpkm[1:10]
		cursor.executemany("update "+tablename+" set "+exp+"=%s where gene = %s",fpkm)
		conn.commit()

def RPKM_DB(cursor,conn,g_t,gene_length,totalreads,samples,intype,insert,tablename):
	sql = 'create table %s select * from %s' %(tablename,gene_length)
	try:
		cursor.execute(sql)
		conn.commit()
		cursor.execute("create index reci on "+tablename+"(gene);")
		conn.commit()
	except:
		print "exists"
	gt = d00.table_2_dict(cursor,g_t,['transc','gene'])
	gl = d00.table_2_dict(cursor,gene_length,['gene','length'])
	for sample in samples:
		total = int(d00.get_sample_info(cursor,sample,totalreads))
		count = {}
		exp = d00.get_sample_info(cursor,sample,'exp')+insert
		try:
			cursor.execute("alter table "+tablename+" add "+exp+" float DEFAULT '0'")
		except:
			print "EXISTS_colume"
		f = open(d00.get_sample_file(cursor,sample,intype))
		for line in f:
			t = re.split('\s+',line)
			if t[0] not in gt:
				continue
			gene = gt[t[0]]
			if gene in gl:
				if gene not in count:
					count[gene] = 0
				count[gene] += int(t[1])
		values = []
		for gene in count:
			values.append([float(count[gene]*1000000*1000)/(total*int(gl[gene])),gene])
		print len(values),values[1]
		cursor.executemany("update "+tablename+" set "+exp+"=%s where gene = %s",values)
		conn.commit()

def COUNT_DB(cursor,conn,g_t,gene_length,samples,intype,insert,tablename):
	sql = 'create table %s select * from %s' %(tablename,gene_length)
	try:
		cursor.execute(sql)
		conn.commit()
		cursor.execute("create index reci on "+tablename+"(gene);")
		conn.commit()
	except:
		print "exists"
	gt = d00.table_2_dict(cursor,g_t,['transc','gene'])
	gl = d00.table_2_dict(cursor,gene_length,['gene','length'])
	for sample in samples:
		count = {}
		exp = d00.get_sample_info(cursor,sample,'exp')+insert
		try:
			cursor.execute("alter table "+tablename+" add "+exp+" float DEFAULT '0'")
		except:
			print "EXISTS_colume"
		f = open(d00.get_sample_file(cursor,sample,intype))
		for line in f:
			t = re.split('\s+',line)
			if t[0] not in gt:
				continue
			gene = gt[t[0]]
			if gene in gl:
				if gene not in count:
					count[gene] = 0
				count[gene] += int(t[1])
		values = []
		for gene in count:
			values.append([count[gene],gene])
		print len(values),values[1]
		cursor.executemany("update "+tablename+" set "+exp+"=%s where gene = %s",values)
		conn.commit()
