import re
import sys
import cPickle as pickle
import subprocess
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
ref = pickle.load(open('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/Ref.all.dat'))
import infra01_pos2info as in1

def blast_fastas(fa_db,fa_target,dbfile,outfile,evalue,wordsize):
	print "HEHE"
	cmd1 = "/data/Analysis/fanxiaoying/software/ncbi-blast-2.2.28+/bin/makeblastdb -dbtype nucl -in %s -out %s " %(fa_db,dbfile)
	cmd2 = "/data/Analysis/fanxiaoying/software/ncbi-blast-2.2.28+/bin/blastn -query %s -task blastn -db %s -outfmt 7 -evalue %s -word_size %s -out %s " %(fa_target,dbfile,evalue,wordsize,outfile)
 	subprocess.call(cmd1,shell=True)
	subprocess.call(cmd2,shell=True)
def blast_fmt7_out_read(file,report):
	file = open(file)
	info = {}
	for line in file:
		if re.match('^#',line):
			continue
		S1 = re.split('\s+',line)
		query_id = S1[0]
		target_id = S1[1]
		match_len = int(S1[3])
		match_percent = round(float(S1[2]))
		query_pos = [int(S1[6]),int(S1[7])]
		target_pos = [int(S1[8]),int(S1[9])]
		if query_id not in info:
			info[query_id] = {}
			info[query_id]['target'] = []
			info[query_id]['match_len'] = []
		if report == 'RC':
			if int(S1[9])<int(S1[8]):
				info[query_id]['target'].append(target_id)
				info[query_id]['match_len'].append(match_len)
		else:
			info[query_id]['target'].append(target_id)
			info[query_id]['match_len'].append(match_len)
	return info
def blast_fmt7_out_read_for_db(file,match_cuf,report):
	file = open(file)
	values = []
	for line in file:
		if re.match('^#',line):
			continue
		S1 = re.split('\s+',line)
		if report == 'RC':
			if int(S1[9])<int(S1[8]) and int(S1[3])>=match_cuf:
				values.append(S1[2:10])
	print values
	return values
def blast_fmt7_out_read_db_miRNA(file,DB_NAME,tablename,report):
	import MySQLdb as mb
	file = open(file)
	values = []
	for line in file:
		if re.match('^#',line):
			continue
		S1 = re.split('\s+',line)
		if report == 'RC':
			if int(S1[9])<int(S1[8]):
				values.append((S1[0],S1[1],S1[6]))
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
	cursor = conn.cursor()
	cursor.executemany("insert into "+tablename+" values(%s,%s,%s) ",values);
	conn.commit()
def blast_genome_positions(species,ranges,evalue,wordsize,report):
	genome = ref[species]['fa']['genome']
	r1_file = './r1.fa'
	r2_file = './r2.fa'
	ranges = in1.get_sequenecs_from_genome_s(species,ranges)
	f = open(r1_file,'w')
	print >>f,">r1\n"+ranges['r1']
	f.close()
	f = open(r2_file,'w')
	print >>f,">r2\n"+ranges['r2']
	f.close()
	blast_fastas(r1_file,r2_file,'./temp_db.db','./temp_blast_r1r2.txt',evalue,wordsize)
	#result = blast_fmt7_out_read_for_db('./temp_blast_r1r2.txt',report)
	#print result
	# if 'r2' not in result:
	# 	out = ['F',0,0]
	# else:
	# 	out = ['T',float(max(result['r2']['match_len']))/(r1[1]-r1[0]),float(max(result['r2']['match_len']))/(r2[1]-r2[0])]
	# return out,result
def blast_genome_multi_positions(species,r1,r2,evalue,wordsize,report):
	genome = ref[species]['fa']['genome']
	r1_file = './r1.fa'
	r2_file = './r2.fa'
	in1.genome_ranges_2_fa_file('mm10',r1,r1_file,'r1')
	in1.genome_ranges_2_fa_file('mm10',r2,r2_file,'r2')
	blast_fastas(r1_file,r2_file,'./temp_db.db','./temp_blast_r1r2.txt',evalue,wordsize)
	result = blast_fmt7_out_read('./temp_blast_r1r2.txt',report)
	return result
def blast_two_sequences(seq1,seq2,evalue,wordsize,report):
	r1_file = open('./r1.fa','w')
	r2_file = open('./r2.fa','w')
	print >>r1_file,'>'+'r1\n'+seq1+'\n'
	print >>r2_file,'>'+'r2\n'+seq2+'\n'
	r1_file.close()
	r2_file.close()
	blast_fastas('./r1.fa','./r2.fa','./temp_db.db','./temp_blast_r1r2.txt',evalue,wordsize)
	result = blast_fmt7_out_read('./temp_blast_r1r2.txt',report)
	return result
def summarize_exons_introns(file,total):
	info = {}
	file = open(file)
	for q_id in range(2,2*total+2):
		query_id = 0.5*q_id
		info[query_id] = {}
		info[query_id]['length'] = 30
		info[query_id]['target'] = []
		info[query_id]['type'] = []
		info[query_id]['qrange'] = []
		info[query_id]['srange'] = []
		info[query_id]['match_len'] = []
		info[query_id]['match_percent'] = []
	print info.keys()
	for line in file:
		if re.match('^#',line):
			continue
		S1 = re.split('\s+',line)
		query_id = float(re.split('#',S1[0])[1])
		target_id = float(re.split('#',S1[1])[1])
		total = int(re.split('#',S1[1])[0])
		match_len = int(S1[3])
		match_percent = round(float(S1[2]))
		query_pos = [int(S1[6]),int(S1[7])]
		target_pos = [int(S1[8]),int(S1[9])]
		info[query_id]['length'] = max(int(S1[3]),info[query_id]['length'])
		info[query_id]['target'].append(target_id)
		info[query_id]['qrange'].append(query_pos)
		info[query_id]['srange'].append(target_pos)
		info[query_id]['match_len'].append(match_len)
		info[query_id]['match_percent'].append(match_percent)
		if target_pos[1]>target_pos[0]:
			info[query_id]['type'].append('Same')
		else:
			info[query_id]['type'].append('RC')
	#print info
	return info