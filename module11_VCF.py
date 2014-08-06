import re, sys, os, copy
import subprocess
import time
import cPickle as pickle
import infra01_pos2info as in1
import MySQLdb as mb
import d00_sample as d00

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
class SNP(dict):
	def __init__(self):
		print "SNP class welcomes you!"
	# def read_VCF_file(self,path,sample_names):
	# 	self['samples'] = sample_names
	# 	file = open(path)
	# 	values = []
	# 	for line in file:
	# 		if re.search('#',line):
	# 			continue
	# 		t = re.split('\s*',line)
	# 		info = re.split(t[7]
	def find_good_quality_SNP_pos(self,group,names,goodsize,QUAL_off,GQ_off,rec):
		self['groupnames'] = names
		self[rec] = {}
		indexs = []
		for i in range(len(names)):
			temp = []
			for j in range(len(group)):
				if group[j] == i:
					temp.append(j)
			indexs.append(temp)	
		for chr in self['chrs']:
			for pos in self[chr]:
				if self[chr][pos]['QUAL'] < QUAL_off:
					continue 
				self[chr][pos]['group_GT'] = ['NA','NA']
				for groupid,i in enumerate(indexs):
					types = []
					number = 0
					for j in i:
						if self[chr][pos]['GQ'][j] >= GQ_off:
							types.append(self[chr][pos]['GT'][j])
					counts = dict([(i, types.count(i)) for i in types])
					GroupType = 'NA'
					for gt in counts:
						if counts[gt] >= goodsize[groupid]:
							GroupType = gt
					self[chr][pos]['group_GT'][groupid] = GroupType
				if 'NA' not in self[chr][pos]['group_GT']:
					counts = dict([(i, types.count(i)) for i in self[chr][pos]['group_GT']])
					if len(counts) == 2:
						if chr not in self[rec]:
							self[rec][chr] = {}
						self[rec][chr][pos] = {}
						self[rec][chr][pos]['GT'] = self[chr][pos]['group_GT']
						self[rec][chr][pos]['ref'] = self[chr][pos]['ref']
						self[rec][chr][pos]['alt'] = self[chr][pos]['alt']
	def get_pos_infos(self,rec,db1,db2):
		poses = copy.deepcopy(self[rec])
		in1.get_infos(db1,db2,poses)
		self[rec] = poses
	def select_target_genes(self,rec,type,genetypes,file):
		outfile = open(file,'w')
		for chr in self[rec]:
			for pos in self[rec][chr]:
				temp = self[rec][chr][pos]
				if self[rec][chr][pos][type]['raw'] == []:
					continue
				if self[rec][chr][pos]['GT'] not in genetypes:
					continue
				print >>outfile,chr,pos,temp['ref'],temp['alt'],temp[type]['genes'][0],temp[type]['transc'][0]
		outfile.close()	
