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
class sampleinfo(dict):
	def __init__(self):
		self['cells'] = []
	def srr_dump_fa(self,ranges,geo,**kwargs):
		self['dump'] = []
		self['merge'] = []
		self['trim'] = []
		if kwargs['pair'] == 'single':
			head = '/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/16.smart_seq/sratool/bin/fastq-dump -gzip'
		else:
			head = '/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/16.smart_seq/sratool/bin/fastq-dump -gzip  -split-3 '
		for i in ranges:
			self['cells'].append(i)
			data_fd = kwargs['data_fd']+'/'+i
			self[i] = {}
			self[i]['files'] = {}
			self[i]['files']['srr'] = []
			self[i]['files']['srr_fq'] = []
			self[i]['files']['raw_fq'] = data_fd+'/'+i+'.fg.gz'
			self[i]['files']['trim_fq'] = data_fd+'/Trim_index.'+i+'.fg.gz'
			self[i]['files']['trim_stats'] = data_fd+'/Trimstat.'+i+'.txt'
			for j in geo[i]['files']:
				SRR = re.split('\.',re.split('/',j)[-1])[0]
				self[i]['files']['srr'].append(re.sub(kwargs['old'],kwargs['data_fd'],j))
				self[i]['files']['srr_fq'].append(data_fd+'/'+SRR+'.fastq.gz')
				if not os.path.exists(data_fd+'/'+SRR+'.fastq.gz'):
					self['dump'].append(head+' -O '+data_fd+' '+re.sub(kwargs['old'],kwargs['data_fd'],j))
			if not os.path.exists(self[i]['files']['raw_fq']):
				cmd = 'cat '
				temp = []
				for j in self[i]['files']['srr_fq']:
					cmd = cmd + ' '+j
					temp.append('rm '+j)
				cmd = cmd+' > '+self[i]['files']['raw_fq']
				self['merge'].append(cmd)
				self['merge'] += temp
			if not os.path.exists(self[i]['files']['trim_fq']):
				pyfile = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/'+kwargs['pyfile']
				cmd = '%s %s %s %s' %(pyfile,self[i]['files']['raw_fq'],self[i]['files']['trim_fq'],self[i]['files']['trim_stats'])
				self['trim'].append(cmd)
	def random_samplings(self,samplename,fq1,fq2,fq3,incount,out_count,datadir,rec):
		if rec not in self:
			self[rec] = []
		script1 = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/random_sample.pair.py '
		script2 = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/random_sample.single.py '
		self['cells'].append(samplename)
		data_dir = datadir+'/'+samplename
		new = samplename
		self[new] = {}
		self[new]['files'] = {}
		self[new]['files']['fq1'] = data_dir+'/SUB_SAMPLE_'+new+'.1.fq.gz' 
		self[new]['files']['fq2'] = data_dir+'/SUB_SAMPLE_'+new+'.2.fq.gz'
		self[new]['files']['fq3'] = data_dir+'/SUB_SAMPLE_'+new+'.S.fq.gz'
		if not os.path.exists(data_dir):
                	os.mkdir(data_dir)
		rate = float(out_count)/incount
		cmd1 = script1+fq1+' '+fq2+' '+self[new]['files']['fq1']+' '+self[new]['files']['fq2']+' '+str(rate)
		cmd2 = script2+fq3+' '+self[new]['files']['fq3']+' '+str(rate)
		if not os.path.exists(self[new]['files']['fq1']):
			self[rec].append(cmd1)
			self[rec].append(cmd2)
	def insert_sample_from_other_db(self,dbname,ids):
		temp = pickle.load(open(dbname))
		for id in ids:
			if id not in self:
				if id not in self['cells']:
					self['cells'].append(id)
				self[id] = temp[id]
	def insert_samples_with_files(self,sample,count,names,pathes):
		self['cells'].append(sample)
		self[sample] = {}
		self[sample]['files'] = {}
		for i in range(0,count):
			self[sample]['files'][names[i]] = pathes[i]
	def process_trim_clean_data(self,A_or_T,samples,outdir,rec):
		self[rec] = []
		if A_or_T == 'A':
			cmd_s = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/p03.trim_A24.py '
		if A_or_T == 'T':
			cmd_s = 'python /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/p03.trim.py '
		for i in samples:
			sample = i
			self[i]['files']['fq1'] = outdir+i+'/TRIM_'+i+'.1.fq.gz'
			self[i]['files']['fq2'] = outdir+i+'/TRIM_'+i+'.2.fq.gz'
			self[i]['files']['fq3'] = outdir+i+'/TRIM_'+i+'.S.fq.gz'
			self[i]['files']['trim_stats'] = outdir+i+'/stat_trim.'+i+'.txt'
			files = self[i]['files']
			cmd = "%s %s %s %s %s %s %s" %(cmd_s,files['raw_fq1'],files['raw_fq2'],files['fq1'],files['fq2'],files['fq3'],files['trim_stats'])
			if not os.path.exists(self[i]['files']['fq1']):
				self[rec].append(cmd)
	def read_trim_info_file(self,samples):
		for i in samples:
			if 'stat' not in self[i]:
				self[i]['stat']={}
			if os.path.exists(self[i]['files']['trim_stats']):
				file = open(self[i]['files']['trim_stats'])
				line = file.readline()
				array = re.split('\s',line)
				self[i]['stat']['total'] = 2*int(array[0])
				self[i]['stat']['trimmed'] = 2*int(array[1])+int(array[2])+int(array[3])
			else:
				print "NO_file",self[i]['files']['trim_stats']
	def process_BWA_pair(self,samples,species,cat,fqs,out,rec):
		self[rec] = []
		sql_command = []
		cmd_scrpit = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BWA.pair.sh'
		for i in samples:
			insert = i+'_'+cat
			outdir = out+'/'+i
			files = self[i]['files']
			files['bam_pair_'+cat] = outdir+'/'+insert+"_pair.bam"
			files['sum_pair_'+cat] = outdir+'/summ_pair.'+insert+'.txt'
			cmd = "%s %s %s %s %s %s %s %s" %(cmd_scrpit,outdir,files[fqs[0]],files[fqs[1]],files['bam_pair_'+cat][:-4],files['sum_pair_'+cat],insert,refall[species]['bwa'][cat])
			if not os.path.exists(self[i]['files']['sum_pair_'+cat]) or os.path.getsize(self[i]['files']['sum_pair_'+cat])<10:
				self[rec].append(cmd)
			sql_command.append([i,'bam_pair_bwa_'+cat,files['bam_pair_'+cat]])
			sql_command.append([i,'sum_pair_'+cat,files['sum_pair_'+cat]])
		return sql_command
	def process_BWA_pair_single(self,samples,species,ref,cat,out,rec):
		self[rec] = []
		sql_command = []
		cmd_s = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BWA.pair_single.sh'
		for i in samples:
			insert = i+cat
			outdir = out+'/'+i
			if not os.path.exists(outdir):
				os.mkdir(outdir)	
			self[i]['files']['bam'+cat] = outdir+'/BAM.'+insert+'.sort.bam'
			self[i]['files']['stat'+cat] = outdir+'/stat.'+insert+'.txt'
			self[i]['files']['sum'+cat] = outdir+'/summ.'+insert+'.txt'
			files = self[i]['files']
			cmd = "%s %s %s %s %s %s %s %s %s %s" %(cmd_s,outdir,files['fq1'],files['fq2'],files['fq3'],files['bam_'+cat][:-4],files['stat'+cat],files['sum'+cat],insert,refall[species]['bwa'][ref])
			if not os.path.exists(self[i]['files']['sum'+cat]):
				self[rec].append(cmd)
			sql_command.append([i,'bam_bwa_'+cat,outdir+'/BAM.'+insert+'.sort.bam'])
			sql_command.append([i,'sum_'+cat,outdir+'/summ.'+insert+'.txt'])
		return sql_command
	def process_TOPCUFF_pair_single(self,samples,species,ref,gtf,cat,out,rec):
		self[rec] = []
		cmd_s = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/TOPCUFF.pair_single.sh '
		for i in samples:
			insert = i+cat
			outdir = out+'/'+i
			if not os.path.exists(outdir):
                                os.mkdir(outdir)
			self[i]['files']['bam'+cat] = outdir+'/BAM.'+insert+'.sort.bam'
                        self[i]['files']['bam_prefix'+cat] = outdir+'/BAM.'+insert+'.sort'
                        self[i]['files']['stat_unmap'+cat] = outdir+'/unmapstat.'+insert+'.txt'
                        self[i]['files']['sum'+cat] = outdir+'/genes.fpkm_tracking'
			files = self[i]['files']
			cmd = "%s %s %s %s %s %s %s %s %s %s" %(cmd_s,outdir,files['fq1'],files['fq2'],files['fq3'],files['bam'+cat],files['stat_unmap'+cat],insert,refall[species]['bowtie2'][ref],refall[species]['gtf'][gtf])
			if not os.path.exists(self[i]['files']['bam'+cat]):
				self[rec].append(cmd)
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
	def read_statfile_refseq(self):
		for i in self['cells']:
			if 'stat' not in self[i]:
				self[i]['stat'] = {}
				print self[i]['files']['stat_refseq']
				file = open(self[i]['files']['stat_refseq'])
				line = file.readline()
				self[i]['stat']['total'] = int(re.split(r'\s',line)[0])
				file.readline()
				line = file.readline()
				self[i]['stat']['mapped'] = int(re.split(r'\s',line)[0])
	def samtools_index(self,bams,rec):
		self[rec] = []
		cmd = 'samtools index '
		for bam in bams:
			for i in self['cells']:
				if not os.path.exists(self[i]['files'][bam]+'.bai'):
					cmd_s = cmd+self[i]['files'][bam]
					self[rec].append(cmd_s)
					print cmd_s

	def SORT_BAM_BY_NAME(self,samples,items,rec):
		self[rec] = []
		for sample in samples:
			for item in items:
				self[sample]['files'][item+'_SORT'] = self[sample]['files'][item][:-4]+'_sort.bam'
				cmd = 'samtools sort -n %s -@ 8 -f %s' %(self[sample]['files'][item],self[sample]['files'][item+'_SORT'])
				self[rec].append(cmd)

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
def UNMAPPED_BWA_Pair_unmapped(cursor,conn,samples,bamtype,folder,rec):
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
def BOWTIE_alignment(cursor,conn,samples,species,ref,ins,outdir,usage,rec):
	cmds = []
	for sample in samples:
		path = outdir+'/BAM.'+sample+'_'+rec+'.bam'
		fq1 = d00.get_sample_file(cursor,sample,ins[0])
		fq2 = d00.get_sample_file(cursor,sample,ins[1])
		cmd = 'bash /data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BOWTIE.pair.sh %s %s %s %s %s' %(fq1,fq2,path[:-4],refall[species]['bowtie2'][ref],usage)
		method = add_files(cmd,"/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/scripts/BOWTIE.pair.sh")
		cursor.execute("replace into files (sample,type,path,method)values(%s,%s,%s,%s)",[sample,rec,path,method])
		conn.commit()
		if not os.path.exists(path):
			cmds.append(cmd)
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
#def BAM_rmd():
#
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
def CUFFLINKS(cursor,conn,species,ref,samples,intype,folder,rec):
	cmds = []
	for sample in samples:
		insert = rec+'_'+sample
		outdir = folder+'/'+insert
		if not os.path.exists(outdir):
			os.mkdir(outdir)
		bam = d00.get_sample_file(cursor,sample,intype)
		path = outdir+'/genes.fpkm_tracking'
		cmd = '/data/Analysis/fanxiaoying/software/cufflinks_221/cufflinks -o %s -p 6 -G %s %s' %(outdir,refall[species]['gtf'][ref],bam)
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
def GENECOUNT_DB(cursor,conn,g_t,gene_length,samples,intype,insert,tablename):
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

def w(a,file):
	f = open(file,'w')
	for i in a:
		print >>f,i
	f.close()

