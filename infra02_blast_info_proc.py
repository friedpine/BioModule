import re,os,sys
import subprocess
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import infra01_pos2info as in01

def blast_fastas(fa_db,fa_target,dbfile,outfile,evalue,wordsize):
	print "Begin Blasting!"
	cmd1 = "/data/Analysis/fanxiaoying/software/ncbi-blast-2.2.28+/bin/makeblastdb -dbtype nucl -in %s -out %s " %(fa_db,dbfile)
 	cmd2 = "/data/Analysis/fanxiaoying/software/ncbi-blast-2.2.28+/bin/blastn -query %s -task blastn -db %s -outfmt 7 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -evalue %s -word_size %s -out %s" %(fa_target,dbfile,evalue,wordsize,outfile)
 	subprocess.call(cmd1,shell=True)
	subprocess.call(cmd2,shell=True)

def blast_fmt7_reading(event,file,match_cuf,report):
	file = open(file)
	values = []
	for line in file:
		if re.match('^#',line):
			continue
		S1 = re.split('\s+',line)
		if report == 'RC':
			if int(S1[9])<int(S1[8]) and int(S1[3])>=match_cuf:
				values.append([event]+S1[0:12])
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


def blast_ref_positions(cursor,species,ref,pos_list1,pos_list2,list1_name,list2_name,evalue,wordsize,match_cuf,folder,record,report,server="TANG"):
	seq1=in01.get_pos_seqs(cursor,species,ref,pos_list1,server)
	seq2=in01.get_pos_seqs(cursor,species,ref,pos_list2,server)
	r1_file = folder+'/'+record+'_r1.fa'
	r2_file = folder+'/'+record+'_r2.fa'
	db_file = folder+'/'+record+'_db.db'
	result_file = folder+'/'+record+'_blast.txt'
	if len(list1_name) != len(seq1) or len(list1_name) != len(seq1):
		print "Names not the Same length with Sequences!"
		return 0
	f = open(r1_file,'w')
	for x in range(len(seq1)):
		f.write('>%s\n%s\n' %(list1_name[x],seq1[x]))
	f.close()
	f = open(r2_file,'w')
	for x in range(len(seq2)):
		f.write('>%s\n%s\n' %(list2_name[x],seq2[x]))
	f.close()
	blast_fastas(r1_file,r2_file,db_file,result_file,evalue,wordsize)
	return blast_fmt7_reading(record,result_file,match_cuf,report)


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
