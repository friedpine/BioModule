import re, sys, os
import subprocess
import time
import cPickle as pickle
import d00_sample as d00
import d02_db_tables as d02
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
refall = pickle.load(open('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/Ref.all.dat'))

def add_files(cmd,file):
	f = open(file)
	out = "#CMD\n"+cmd+"\n"+"#FILE\t"+file+"\n"
	for i in f:
		out+=i
	return out

def update_STATA(cursor,conn):
	cursor.execute("select sample,type,path from files")
	for r0 in cursor.fetchall():
		try:
			s = os.path.getsize(r0[2])
			cursor.execute("update files set state=%s where sample=%s and type=%s",(s,r0[0],r0[1]))
		except:
			print r0
	conn.commit()

def circularRNA_tophat_pair(cursor,conn,samples,species,ref,out,in1,in2,rec):
	cmds = []
	for sample in samples:
		outdir = out+'/'+rec+"_"+sample
		if not os.path.exists(outdir):
			os.mkdir(outdir)
		path = outdir+'/'+rec+"_"+sample+'.bam'
		c = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/cRNA_TOPHAT.pair.sh '
		fq1 = d00.get_sample_file(cursor,sample,in1)
		fq2 = d00.get_sample_file(cursor,sample,in2)
		cmd = "%s %s %s %s %s %s %s" %(c,outdir,fq1,fq2,path[:-4],rec,refall[species]['bowtie2'][ref])
		if not os.path.exists(path):
			cmds.append(cmd)
		cursor.execute("insert ignore into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,cmd])
		conn.commit()
	return cmds

def COMMAND_generator(cursor,conn,samples,template,infiles,folder,suffix,rec):
	cmds = []
	for sample in samples:
		lists = []
		for i in infiles:
			lists.append(d00.get_sample_file(cursor,sample,i))
		path = folder+'/'+rec+'_'+sample+suffix
		if suffix=='/':
			path = path[:-1]
		if suffix=='/' and not os.path.exists(path):
			os.mkdir(path)
		lists.append(path)
		cmd = template
		for i in range(len(lists)):
			cmd = cmd.replace("#"+str(i),lists[i])
		cmd = cmd.replace("#sample",sample)
		if not (rec == '' or rec == 'NA'):
			cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,cmd])
		cmds.append(cmd)
	conn.commit()
	return cmds

def UNMAPPED_Pair_unmapped(cursor,conn,samples,bamtype,folder,rec):
	cmds = []
	for sample in samples:
		bam = d00.get_sample_file(cursor,sample,bamtype)
		path = folder+'/'+rec+'_'+sample+'.fa.gz'
		cmd = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/assembly/unmapped_reads_BWA.py ' +bam+' '+path+' '+sample
		method = add_files(cmd,"/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/assembly/unmapped_reads_BWA.py")
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,method])
		cmds.append(cmd)
	conn.commit()
	return cmds

def FILES_GROUPER(cursor,conn,samples,newsName,intype,sep):
	newname = ""
	for sample in samples:
		newname += d00.get_sample_file(cursor,sample,intype)+sep
	print newname[:-1]
	method = "FILES_GROUPER #"+" "+intype
	cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[newsName,intype,newname,method])
	conn.commit()

def MAPPED_SINGLE(cursor,conn,samples,bamtype,folder,rec):
	cmds = []
	for sample in samples:
		bam = d00.get_sample_file(cursor,sample,bamtype)
		path = folder+'/'+rec+'_'+sample+'.fa.gz'
		cmd = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/assembly/mapped_reads_single.py ' +bam+' '+path+' '+sample
		method = add_files(cmd,"/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/assembly/mapped_reads_single.py")
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,method])
		cmds.append(cmd)
	conn.commit()
	return cmds

def BOWTIE_PAIRED(cursor,conn,samples,species,ref,ins,outdir,usage,rec):
	cmds = []
	for sample in samples:
		path = outdir+'/BAM.'+rec+'_'+sample+'.bam'
		fq1 = d00.get_sample_file(cursor,sample,ins[0])
		fq2 = d00.get_sample_file(cursor,sample,ins[1])
		cmd = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BOWTIE.pair.sh %s %s %s %s %s' %(fq1,fq2,path[:-4],refall[species]['bowtie2'][ref],usage)
		method = add_files(cmd,"/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BOWTIE.pair.sh")
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,method])
		conn.commit()
		if not os.path.exists(path):
			cmds.append(cmd)
	return cmds

def BOWTIE_SINGLE(cursor,conn,samples,species,ref,ins,outdir,usage,rec):
	cmds = []
	return cmds

def BWA_SINGLE(cursor,conn,specise,ref,samples,intype,folder,rec):
	cmds = []
	for sample in samples:
		insert = rec+'_'+sample
		fq = d00.get_sample_file(cursor,sample,intype)
		path = folder+'/BAM'+insert+'.bam'
		cmd = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BWA.single.sh ' +folder+' '+fq\
				+' '+refall[specise]['bwa'][ref]+" "+path[:-4]+" "+insert
		method = add_files(cmd,"/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BWA.single.sh")
		cursor.execute("insert ignore into files (sample,type,path,method)values(%s,%s,%s,%s) ",[sample,rec,path,method])
		cmds.append(cmd)
	conn.commit()
	return cmds

def BWA_PAIRED(cursor,conn,specise,ref,samples,intype,folder,rec):
	cmds = []
	for sample in samples:
		insert = rec+'_'+sample
		fq1 = d00.get_sample_file(cursor,sample,intype[0])
		fq2 = d00.get_sample_file(cursor,sample,intype[1])
		path = folder+'/BAM'+insert+'.bam'
		cmd = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BWA.pair.sh ' +folder+' '+fq1+' '+fq2\
				+' '+refall[specise]['bwa'][ref]+" "+path[:-4]+" "+insert
		method = add_files(cmd,"/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BWA.pair.sh")
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s) ",[sample,rec,path,method])
		cmds.append(cmd)
	conn.commit()
	return cmds

def TOPHAT_SINGLE(cursor,conn,specise,ref,samples,intype,report,folder,rec):
	cmds = []
	for sample in samples:
		insert = rec+'_'+sample
		outdir = folder+'/'+insert
		bam = d00.get_sample_file(cursor,sample,intype)
		path = outdir+'/BAM'+insert+'.bam'
		cmd = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/TOPHAT.single.sh ' +outdir+' '+bam\
				+' '+refall[specise]['bowtie2'][ref]+" "+path[:-4]+" "+insert+" "+report
		method = add_files(cmd,"/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/TOPHAT.single.sh")
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,method])
		cmds.append(cmd)
	conn.commit()
	return cmds

def TOPHAT_PAIRED(cursor,conn,specise,ref,samples,intype,report,folder,rec):
	cmds = []
	for sample in samples:
		insert = rec+'_'+sample
		outdir = folder+'/'+insert
		fq1 = d00.get_sample_file(cursor,sample,intype[0])
		fq2 = d00.get_sample_file(cursor,sample,intype[1])
		path = outdir+'/BAM'+insert+'.bam'
		cmd = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/TOPHAT.pair.sh ' +outdir+' '+fq1+' '+fq2\
				+' '+refall[specise]['bowtie2'][ref]+" "+path[:-4]+" "+insert+" "+report
		method = add_files(cmd,"/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/TOPHAT.pair.sh")
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,method])
		cmds.append(cmd)
	conn.commit()
	return cmds

def BAM_FLAGSTAT(cursor,conn,samples,intype,folder,rec):
	cmds = []
	for sample in samples:
		path = folder+'/flagstat.'+sample+"_"+intype+'.txt'
		f1 = d00.get_sample_file(cursor,sample,intype)
		cmd = 'samtools flagstat '+f1+' > '+path
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,cmd])
		cmds.append(cmd)
	conn.commit()
	return cmds

def Read_BWA_flagstat(cursor,conn,samples,intype,recs):
	d02.check_table_colume(cursor,conn,'samples',recs[0],'INT')
	d02.check_table_colume(cursor,conn,'samples',recs[1],'INT')
	for sample in samples:
		f1 = open(d00.get_sample_file(cursor,sample,intype))
		a = re.split('\s+',f1.readline())[0]
		f1.readline()
		b = re.split('\s+',f1.readline())[0]
		sql = "update %s set %s=%s,%s=%s where sample='%s'" %('samples',recs[0],a,recs[1],b,sample)
		print sql
		cursor.execute(sql)
	conn.commit()

def SUMMARIZE(cursor,conn,samples,intype,folder,para,rec):
	cmds = []
	for sample in samples:
		insert = rec+'_'+sample
		bam = d00.get_sample_file(cursor,sample,intype)
		path = folder+'/summarize.'+insert+'.txt'
		cmd = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/SUMMARIZE.sh '+bam+' '+path+' '+para
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,cmd])
		cmds.append(cmd)
	conn.commit()
	return cmds


def w(a,file):
	f = open(file,'w')
	for i in a:
		print >>f,i
	f.close()

