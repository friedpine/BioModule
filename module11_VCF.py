import re, sys, os, copy
import subprocess
import time
import cPickle as pickle
import infra01_pos2info as in1
import MySQLdb as mb
import d00_sample as d00
import module00_samples as m00

def TEST():
	print "TEST"

def read_VCF_file(cursor1,SAMPLE_IN,type,DB_NAME,tablename,limit,counts,samples):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
	cursor = conn.cursor()
	sample_infos = ''
	for sample in samples:
		sample_info = " %s_G varchar(5) DEFAULT NULL,%s_0 int(11) DEFAULT NULL,%s_1 int(11) DEFAULT NULL,%s_2 int(11) DEFAULT NULL,%s_DP int(11) DEFAULT NULL,%s_GQ int(11) DEFAULT NULL," %(sample,sample,sample,sample,sample,sample)
		sample_infos += sample_info
	sql = """CREATE TABLE %s (
	  `chr` varchar(20) NOT NULL DEFAULT '',
	  `pos` int(11) NOT NULL DEFAULT '0',
	  `Ref` varchar(30) DEFAULT NULL,
	  `Alt` varchar(30) NOT NULL DEFAULT '',
	  `Qual` float DEFAULT NULL,
	  `DP` int(11) DEFAULT NULL,
	  `DP_alt` int(11) DEFAULT NULL,
	  `FQ` float DEFAULT NULL,
	  `AF1` float DEFAULT NULL,
	  `AC1` float DEFAULT NULL,
	  %s
	  PRIMARY KEY (`chr`,`pos`,`Alt`)
	) ENGINE=InnoDB DEFAULT CHARSET=latin1""" %(tablename,sample_infos)
	try:
		cursor.execute(sql)
	except:
		print "EXISTS"
	print cursor1,sample,type
	path=d00.get_sample_file(cursor1,SAMPLE_IN,type)
	file = open(path)
	values = []
	for line in file:
		if re.search('#',line):
			continue
		t = re.split('\s*',line)
		info = {} 
		for i in re.split(';',t[7]):
			a = re.split('=',i)
			if len(a)>1:
				info[a[0]] = a[1]
		if len(t[3])>limit:
			t[3]=t[3][0:20]
			continue
		if len(t[4])>limit:
			t[4]=t[4][0:limit]
			continue
		if "DP4" in info:
			DP4 = sum([int(i) for i in re.split(",",info['DP4'])][2:4])		
		else:
			DP4=0
		value = (t[0],t[1],t[3],t[4],t[5],info['DP'],DP4,info['FQ'],info['AF1'],info['AC1'])
		for i in range(counts):
			value += tuple(re.split(':|,',t[9+i]))
		if len(value)!=10+counts*6:
			a = 10
		else:
			values.append(value)
	cmd = "insert into "+tablename+" values(%s"+",%s"*(9+counts*6)+")"
	print len(values)
	cursor.executemany(cmd,values);
	conn.commit()
	cursor.close()
	conn.close()
def read_VCF_file_single(cursor,conn,DB_NAME,tablename,samples,type):
	limit = 30
	sample_infos = ''
	for sample in samples:
		sample_info = " %s_DP varchar(5) DEFAULT '0',%s_alt float DEFAULT '0'," %(sample,sample)
		sample_infos += sample_info
	sql = """CREATE TABLE %s (
	  `chr` varchar(20) NOT NULL DEFAULT '',
	  `pos` int(11) NOT NULL DEFAULT '0',
	  `Ref` varchar(30) DEFAULT NULL,
	  `Alt` varchar(30) NOT NULL DEFAULT '',
	  %s
	  PRIMARY KEY (`chr`,`pos`,`Alt`)
	) ENGINE=InnoDB DEFAULT CHARSET=latin1""" %(tablename,sample_infos)
	print sql
	try:
		cursor.execute(sql)
	except:
		print "EXISTS"
	for sample in samples:
		path = d00.get_sample_file(cursor,sample,type)
		file = open(path)
		values = []
		for line in file:
			if re.search('#',line):
				continue
			t = re.split('\s*',line)
			info = {} 
			for i in re.split(';',t[7]):
				a = re.split('=',i)
				if len(a)>1:
					info[a[0]] = a[1]
			if 'DP4' not in info:
				continue
			DP4 = re.split(',',info['DP4'])
			if len(t[3])>limit:
				t[3]=t[3][0:limit]
				continue
			if len(t[4])>limit:
				t[4]=t[4][0:limit]
				continue
			value = (t[0],t[1],t[3],t[4],info['DP'],float(int(DP4[2])+int(DP4[3]))/int(info['DP']))
			values.append(value)
		cmd = "insert into %s (chr,pos,Ref,Alt,%s,%s)values(%%s,%%s,%%s,%%s,%%s,%%s) on duplicate key update %s=values(%s),%s=values(%s)" %(tablename,sample+'_DP',sample+'_alt',sample+'_DP',sample+'_DP',sample+'_alt',sample+'_alt')
		print cmd,values[0]
		cursor.executemany(cmd,values)
		conn.commit()
	cursor.close()
	conn.close()
def CELLS_GROUPS_VCF(cursor,conn,samples,group_name,file_type,outdir,cmd_file):
	m00.FILES_GROUPER(cursor,conn,samples,group_name,file_type,' ')
	w1 = m00.COMMAND_generator(cursor,conn,[group_name],"samtools mpileup -Duf /data/Analysis/fanxiaoying/database/hg19/00.genome/genome.fa #0 | /data/Analysis/fanxiaoying/software/samtools-0.1.19/bcftools/bcftools view -bvcg -> #1",[file_type],outdir,'.bcf','SNV')
	w2 = m00.COMMAND_generator(cursor,conn,[group_name],"/data/Analysis/fanxiaoying/software/samtools-0.1.19/bcftools/bcftools view #0 | perl /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/vcfutils.pl varFilter -d 5 >#1",['SNV'],outdir,'.vcf','SNV_d5')
	m00.w(w1+w2,cmd_file)
	
def Single_Sample_SNV_FILE(cursor,conn,sample,fq_type,bam_type,bam_type_rmd,SNV_type,SNV_process_type,outdir,cmd_file):
	map_cmd = m00.BWA_PAIRED(cursor,conn,"hg19","genome",[sample],fq_type,outdir,bam_type)
	rmdup = m00.COMMAND_generator(cursor,conn,[sample],"samtools rmdup -S #0 #1",[bam_type],outdir,'.bam',bam_type_rmd)
	mpileup = m00.COMMAND_generator(cursor,conn,[sample],"samtools mpileup -Duf /data/Analysis/fanxiaoying/database/hg19/00.genome/genome.fa #0 | /data/Analysis/fanxiaoying/software/samtools-0.1.19/bcftools/bcftools view -bvcg -> #1",[bam_type_rmd],outdir,'.bcf',SNV_type)
	process = m00.COMMAND_generator(cursor,conn,[sample],"/data/Analysis/fanxiaoying/software/samtools-0.1.19/bcftools/bcftools view #0 | perl /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/vcfutils.pl varFilter -d 5 >#1",[SNV_type],outdir,'.vcf',SNV_process_type)
	m00.w(map_cmd+rmdup+mpileup+process,cmd_file)
	
