from __future__ import division
import sys,os,re,string,time
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
import pairwise2
from pairwise2 import format_alignment
import subprocess
import d00_sample as d00
import cPickle as pickle
import gzip
import numpy as np
import MySQLdb as mb
import matplotlib
#matplotlib.use('Cairo')
import matplotlib.pyplot as plt
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
Ref = pickle.load(open('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/Ref.all.dat'))

def complementary(seq):
	rc = ''
	for i in seq:
		if i == 'A':
			rc=rc+'T'
		elif i == 'T':
			rc=rc+'A'
		elif i == 'C':
			rc=rc+'G'
		elif i == 'G':
			rc=rc+'C'
		else:
			rc=rc+i
	return rc
def makes_anchors_fq(cursor,conn,samples,bamfile,outdir,length,insert,rec):
	sql = []
	cmds = []
	for sample in samples:
		anchor1 = outdir+'/Anchor.'+insert+'_'+sample+'.1.fq.gz'
		anchor2 = outdir+'/Anchor.'+insert+'_'+sample+'.2.fq.gz'
		if os.path.exists(anchor1):
			print "ALREADY EXISTS"
		else:
			bam_pe = d00.get_sample_file(cursor,sample,bamfile) 
			cmd = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/cRNA/make_anchors.py %s %s %s %s' %(bam_pe,anchor1,anchor2,length)
			cmds.append(cmd)
		sql.append([sample,rec[0],anchor1,cmd])
		sql.append([sample,rec[1],anchor2,cmd])
	cursor.executemany("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",sql)
	conn.commit()
	return cmds
def BOWTIE_alignment(cursor,conn,samples,species,in1,in2,rec,para,outdir):
	cmds = []
	for sample in samples:
		outsam = outdir+'/SAM.anchor.'+sample+'_'+rec+'.sam'
		path = outdir+'/BAM.anchor_genome.'+sample+'_'+rec+'.bam'
		fq1 = d00.get_sample_file(cursor,sample,in1)
		fq2 = d00.get_sample_file(cursor,sample,in2)
		cmd = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/cRNA/BOWTIE.pair.sh %s %s %s %s %s %s' %(fq1,fq2,outsam,path[:-4],Ref[species]['bowtie2']['genome'],para)
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,cmd])
		conn.commit()
		if not os.path.exists(path):
			cmds.append(cmd)
	return cmds
def summarize_pair_map_conditions_db(cursor,samples,dbname,tablename,in1,dist,anchor_length):
	cmds = []
	for sample in samples:
		bamfile = d00.get_sample_file(cursor,sample,in1)
		cmd = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/cRNA/summarize_reads_db.py %s %s %s %s %s %s' %(dbname,tablename,bamfile,dist,sample,anchor_length)
		cmds.append(cmd)
	return cmds
def candicate_reads_and_mate_map_2_genome_refseq(cursor,conn,sample,species,reftype,folder_raw,db_name,tablename):
	sfolder = folder_raw+'/circ_SELF_'+sample
	mfolder = folder_raw+'/circ_MATE_'+sample
	if not os.path.exists(sfolder):
		os.mkdir(sfolder)
	if not os.path.exists(mfolder):
		os.mkdir(mfolder)
	sams = sfolder+'/temp_bowtie.sam'
	samm = mfolder+'/temp_bowtie.sam'
	mate_fa = mfolder+'/Mate_seq.fa'
	mate_genome = mfolder+'/Mate_bowtie_'+species+'.bam'
	mate_refseq = mfolder+'/Mate_bowtie_refseq.bam'
	self_fa = sfolder+'/Self_seq.fa'
	self_genome = sfolder+'/Self_bowtie_'+species+'.bam'
	self_refseq = sfolder+'/Self_bowtie_refseq.bam'
	if not os.path.exists(mfolder+'/Mate_seq.fa') or os.path.getsize(mfolder+'/Mate_seq.fa')==0:
		f1 = open(mate_fa,'w')
		f2 = open(self_fa,'w')
		conn1=mb.connect(host="localhost",user="root",passwd="123456",db=db_name)
		cursor1 = conn1.cursor()
		cursor1.execute("select id,seq,mate from "+tablename)
		results = cursor1.fetchall()
		for t in results:
			f2.write('>'+str(t[0])+"\n"+t[1]+"\n")
			f1.write('>'+str(t[0])+"\n"+t[2]+"\n")
		f2.close()
	values = [[sample,'mate_genome',mate_genome],[sample,'mate_refseq',mate_refseq],[sample,'self_genome',self_genome],[sample,'self_refseq',self_refseq]]
	cursor.executemany("insert ignore into files values(%s,%s,%s,NULL)",values)
	conn.commit()
	cmd_bwa = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BWA.single_simp.sh '
	cmd_top = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/cRNA_bowtie_single_refseq_genome.sh '
	cmd1 = cmd_bwa+Ref[species]['bwa'][reftype]+' '+self_fa+' '+sfolder+' '+self_refseq[:-4]+' self_'+reftype
	cmd2 = cmd_bwa+Ref[species]['bwa'][reftype]+' '+mate_fa+' '+mfolder+' '+mate_refseq[:-4]+' mate_'+reftype
	cmd3 = cmd_top+Ref[species]['bowtie2']['genome']+' '+self_fa+' '+sfolder+' '+self_genome[:-4]+' '+Ref[species]['gtf']['ensembl']
	cmd4 = cmd_top+Ref[species]['bowtie2']['genome']+' '+mate_fa+' '+mfolder+' '+mate_genome[:-4]+' '+Ref[species]['gtf']['ensembl']
	print cmd1,'\n',cmd2,'\n',cmd3,'\n',cmd4
def read_self_mate_map_data_and_trim(cursor,conn,sample,species,DB_NAME,tablename,table_raw):
	conn1=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
	cursor1 = conn1.cursor()
	files = ['self_refseq','mate_refseq']
	rec = ['S_R','M_R']
	for index,filename in enumerate(files):
		values = []
		cmd = 'samtools view -F 4 '+d00.get_sample_file(cursor,sample,filename)
		p1 = subprocess.Popen(cmd,shell = True,stdout=subprocess.PIPE)
		for line in p1.stdout:
			t = str.split(line)
			values.append((t[0],rec[index],'T','101M'))
		cursor1.executemany("insert into "+tablename+" values(%s,%s,%s,1000,%s) ",values);
		conn1.commit()
	files = ['mate_genome']
	rec = ['M_G']
	for index,filename in enumerate(files):
		values = []
		cmd = 'samtools view -F 4 '+d00.get_sample_file(cursor,sample,filename)
		p1 = subprocess.Popen(cmd,shell = True,stdout=subprocess.PIPE)
		for line in p1.stdout:
			t = str.split(line)
			if species == 'mm10':
				t[2] = 'chr'+t[2]
			if len(t[5])>20:
				t[5] = t[5][0:20]
			values.append((t[0],rec[index],t[2],t[3],t[5]))
		cursor1.executemany("insert into "+tablename+" values(%s,%s,%s,%s,%s) ",values);
		conn1.commit()
def summarize_realignmap_data(DB_NAME,table_raw,table_realign):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
	cursor = conn.cursor()
	cursor.execute("select id from "+table_realign+" where rec= 'M_R'")
	results = cursor.fetchall()
	values = []
	for i in results:
		values.append((1,i[0]))
	cursor.executemany("update "+table_raw+" set mate_gencode= %s where id = %s ",values);
	conn.commit()
	cursor.execute("select id from "+table_realign+" where rec= 'S_R'")
	results = cursor.fetchall()
	values = []
	for i in results:
		values.append((1,i[0]))
	cursor.executemany("update "+table_raw+" set self_gencode= %s where id = %s ",values);
	conn.commit()
	cursor.close()
	conn.close()

class cirrna_meth1(dict):
	def __init__(self,sample,info):
		self['event'] = {}
		self['event_info'] = {}
		self['sample'] = sample
		self['info'] = info
	def makes_anchors_fq(self,samples,outdir,length,rec):
		self[rec] = []
		for sample in samples:
			files = self['info'][sample]['files']
			anchor1 = outdir+'/Anchor.'+sample+'.1.fq.gz'
			anchor2 = outdir+'/Anchor.'+sample+'.2.fq.gz'
			files['anchor1'] = anchor1
			files['anchor2'] = anchor2
			if os.path.exists(anchor1):
				print "ALREADY EXISTS"
			else:
				bam_pe = files['BAM_PE_cRNA']
				cmd = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/cRNA_making_anchors_pair.py %s %s %s %s' %(bam_pe,anchor1,anchor2,length)
				self[rec].append(cmd)
	def BOWTIE_alignment(self,samples,species,outdir,rec):
		self[rec] = []
		for sample in samples:
			outsam = outdir+'/SAM.anchor.'+sample+'.sam'
			bam_prefix = outdir+'/BAM.anchor.'+sample
			files = self['info'][sample]['files']
			files['ANCHOR_BAM_SORT'] = bam_prefix+'.bam'
			cmd = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BOWTIE.pair.sh %s %s %s %s %s' %(files['anchor1'],files['anchor2'],outsam,bam_prefix,Ref[species]['bowtie2']['genome'])
			if not os.path.exists(bam_prefix+'.bam'):
				self[rec].append(cmd)	
	def summarize_pair_map_conditions_db(self,samples,dbname,tablename,dist,record,anchor_length):
		self[record] = []
		for sample in samples:
			files = self['info'][sample]['files']
			cmd = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/cRNA_summarize_reads_db.py %s %s %s %s %s' %(dbname,tablename,files['ANCHOR_BAM_SORT'],dist,sample,anchor_length)
			self[record].append(cmd)
	def candicate_reads_and_mate_map_2_genome_refseq(self,samples,species,folder_raw,db_name):
		if 'files' not in self:
			self['files'] = {}
		sfolder = folder_raw+'/SELF_info'
		mfolder = folder_raw+'/MATE_info'
		if not os.path.exists(sfolder):
			os.mkdir(sfolder)
		if not os.path.exists(mfolder):
			os.mkdir(mfolder)
		sams = sfolder+'/temp_bowtie.sam'
		samm = mfolder+'/temp_bowtie.sam'
		self['files']['mate_fa'] = mfolder+'/Mate_seq.fa'
		self['files']['mate_genome'] = mfolder+'/Mate_bowtie_'+species+'.bam'
		self['files']['mate_refseq'] = mfolder+'/Mate_bowtie_refseq.bam'
		self['files']['self_fa'] = sfolder+'/Self_seq.fa'
		self['files']['self_genome'] = sfolder+'/Self_bowtie_'+species+'.bam'
		self['files']['self_refseq'] = sfolder+'/Self_bowtie_refseq.bam'
		if not os.path.exists(mfolder+'/Mate_seq.fa'):
			f1 = open(self['files']['mate_fa'],'w')
			f2 = open(self['files']['self_fa'],'w')
			conn=mb.connect(host="localhost",user="root",passwd="123456",db=db_name)
			cursor = conn.cursor()
			cursor.execute("select id,seq,mate from raw1")
			results = cursor.fetchall()
			for t in results:
				f2.write('>'+str(t[0])+"\n"+t[1]+"\n")
				f1.write('>'+str(t[0])+"\n"+t[2]+"\n")
			f2.close()
		cmd_bwa = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BWA.single_simp.sh '
		cmd_top = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/cRNA_bowtie_single_refseq_genome.sh '
		cmd1 = cmd_bwa+Ref[species]['bwa']['gencode']+' '+self['files']['self_fa']+' '+sfolder+' '+self['files']['self_refseq'][:-4]+' self_gencode'
		cmd2 = cmd_bwa+Ref[species]['bwa']['gencode']+' '+self['files']['mate_fa']+' '+mfolder+' '+self['files']['mate_refseq'][:-4]+' mate_gencode'
		cmd3 = cmd_top+Ref[species]['bowtie2']['genome']+' '+self['files']['self_fa']+' '+sfolder+' '+self['files']['self_genome'][:-4]+' '+Ref[species]['gtf']['ensembl']
		cmd4 = cmd_top+Ref[species]['bowtie2']['genome']+' '+self['files']['mate_fa']+' '+mfolder+' '+self['files']['mate_genome'][:-4]+' '+Ref[species]['gtf']['ensembl']
		print cmd1,'\n',cmd2,'\n',cmd3,'\n',cmd4
	def read_self_mate_map_data_and_trim(self,species,DB_NAME):
		files = ['self_refseq','mate_refseq']
		rec = ['S_R','M_R']
		# for index,filename in enumerate(files):
		# 	values = []
		# 	cmd = 'samtools view -F 4 '+self['files'][filename]
		# 	p1 = subprocess.Popen(cmd,shell = True,stdout=subprocess.PIPE)
		# 	for line in p1.stdout:
		# 		t = str.split(line)
		# 		values.append((t[0],rec[index],'T','101M'))
		# 	conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
		# 	cursor = conn.cursor()
		# 	cursor.executemany("""insert into realign values(%s,%s,%s,1000,%s) """,values);
		# 	conn.commit()
		# files = ['self_genome','mate_genome']
		# rec = ['S_G','M_G']
		# for index,filename in enumerate(files):
		# 	values = []
		# 	cmd = 'samtools view -F 4 '+self['files'][filename]
		# 	p1 = subprocess.Popen(cmd,shell = True,stdout=subprocess.PIPE)
		# 	for line in p1.stdout:
		# 		t = str.split(line)
		# 		if species == 'mm10':
		# 			t[2] = 'chr'+t[2]
		# 		if len(t[5])>20:
		# 			t[5] = t[5][0:20]
		# 		values.append((t[0],rec[index],t[2],t[3],t[5]))
		# 	conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
		# 	cursor = conn.cursor()
		# 	cursor.executemany("""insert into realign values(%s,%s,%s,%s,%s) """,values);
		# 	conn.commit()
		conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
		cursor = conn.cursor()
		print "GET_IDs",time.time()
		cursor.execute("select id from realign where rec= 'M_R'")
		results = cursor.fetchall()
		values = []
		print "GOT_IDs",time.time()
		for i in results:
			values.append((1,i[0]))
		print "values is ready",len(values),time.time()
		cursor.executemany("""update raw1 set mate_gencode= %s where id = %s """,values);
		print "EXECUTE_DOWN",time.time()
		conn.commit()
		cursor.close()
		conn.close()
	def connect_anchors_with_konwn_exons(self,samples,folder,exondb,map2genome_up,map2refseq_up,precision,rec,species):
		print "HEHE"
	
#	def events_full_length_info(self)	
	def run_pipeline(self,cmdname,n):
		handle = []
		for i in range(0,n):
			handle.append('PP'+str(i))
		info =  self[cmdname]
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
					break	
def integrate_events(DB_NAME,table1,table2):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
	cursor = conn.cursor()
	cursor.execute("select * from "+table1+" where total_ex<55 and total_ex >45")
	r = cursor.fetchall()
	values = []
	for i in r:
		values.append((i[2],i[3],i[4],i[5]))
	cursor.executemany("insert ignore into "+table2+"(id,transc,l_ex,r_ex)values(%s,%s,%s,%s)",values)
	conn.commit()
def events_full_length_info(DB_NAME,db_ref,table1):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
	cursor = conn.cursor()
	cursor.execute("select * from "+table1)
	values = []
	r0 = cursor.fetchall()
	for i in r0:
		cursor.execute("select * from "+db_ref+".t_ex where transc = %s and ex_order >= %s and ex_order<= %s",[i[1],i[2],i[3]])
		r1 = cursor.fetchall()
		length = 0
		sites = 0
		for j in r1:
			length += j[5]-j[4]
			cursor.execute("select miRNA_count from "+db_ref+".exons where chr=%s and left_pos = %s and right_pos = %s",[j[2],j[4],j[5]])
			sites += int(cursor.fetchall()[0][0])
		values.append((length,sites,i[0]))
	cursor.executemany("update "+table1+" set length = %s,miRNA_sites = %s where id = %s",values)
	conn.commit()
def events_type(DB_NAME,db_ref,table1):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
	cursor = conn.cursor()
	cursor.execute("select * from "+table1)
	r0 = cursor.fetchall()
	values = []
	for i in r0:
		cursor.execute("select cate from "+db_ref+".g_t join "+db_ref+".gene where transc = %s and g_t.ENSG = gene.ENSG",[i[1]])
		r1 = cursor.fetchall()
		values.append((r1[0][0],i[0]))
	cursor.executemany("update "+table1+" set gene_type = %s where id = %s",values)
	conn.commit()
def events_expression_level(DB_NAME,table1,raw1,samples):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
	cursor = conn.cursor()
	for sample in samples:
		sample = re.sub('-','_',sample) 
		try:
			cursor.execute("alter table "+table1+" add "+'mm_'+str(sample)+" INT default 0")
			conn.commit()
		except:
		 	print "ALREADY EXISTS"
	cursor.execute("select id from all_events")
	r0 = cursor.fetchall()
	events = []
	for i in r0:
		print i[0]
		cursor.execute("insert into all_events_expression (id) values (%s)",[i[0]])
	conn.commit()
	for sample in samples:
		sample = re.sub('-','_',sample) 
		cursor.execute("select event,count(*) from "+raw1+" where sample = %s group by event",[sample])
		r1 = cursor.fetchall()
		values = [(i[1],i[0]) for i in r1]
		print values
		cursor.executemany("update "+table1+" set "+'mm_'+str(sample)+" = %s where id = %s",values)
		conn.commit()
#def events_expression_level_genes()	

	
