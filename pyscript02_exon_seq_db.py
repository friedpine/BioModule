import os,subprocess,sys
import module10_make_database as m10
import cPickle as pickle
import time,string,re
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
Ref = pickle.load(open('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/Ref.all.dat'))

##CONFIGURATIONS
spec = 'mm10'
regions = ['exon','intron']
datadb = "Ref.intron_exon.mm10.dat"
gtf_type = 'refseq'

exondb = m10.exon_database()
if 'exon' in regions:
	bedfile = Ref[spec]['bed']['exons']
	exon_fa = './temp_'+spec+'_exon.fa'
	exondb.gtf_exon2bed(Ref[spec]['gtf'][gtf_type],bedfile)
	cmd_bed = "/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/fastaFromBed -fi %s -bed %s -fo %s" %(Ref[spec]['fa']['genome'],Ref[spec]['bed']['exons'],exon_fa)
	if not os.path.exists(exon_fa):
		subprocess.call(cmd_bed,shell=True)
	exondb.read_exon_fasta_file(exon_fa)
	exondb.insert_transcript_exon_count_info('exon_info')
	pickle.dump(exondb,open(datadb,"wb"),True)
if 'intron' in regions:
	bedfile_50 = Ref[spec]['bed']['introns_50']
	bedfile_all = Ref[spec]['bed']['introns']
	exondb.gtf_intron2bed(Ref[spec]['gtf'][gtf_type],bedfile_50,bedfile_all,49)
	cmd_bed = "/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/fastaFromBed -fi %s -bed %s -fo %s" %(genome,bedfile_50,intron_fa)
	if not os.path.exists(intron_fa):
		subprocess.call(cmd_bed,shell=True)
	exondb.read_intron_fasta_file('intron_info',intron_fa)
	pickle.dump(exondb,open(datadb,"wb"),True)
