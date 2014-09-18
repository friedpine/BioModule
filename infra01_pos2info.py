import re,copy,sys
import cPickle as pickle
import subprocess
import infra00_ranges_operate as in0
import mmap
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
#ref = pickle.load(open('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/Ref.all.dat'))

def reverse_complementary(seq):
	rc = ''
	reverse = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	for i in seq:
		rc=rc+reverse[i]
	return rc[::-1]


def get_seq_mm_mmap(chr,start,end):
	f = mmap_fasta("/data/Analysis/fanxiaoying/database/mm10/01.bowtie/"+chr+'.fa')
	return f.getseq(start-1,end)

def get_seq_hg_mmap(chr,start,end):
	f = mmap_fasta("/data/Analysis/fanxiaoying/database/hg19/00.genome/"+chr+'.fa')
	return f.getseq(start-1,end).upper()

def get_multiseqs_mmap(species,ranges):
	outs = ['N']*len(ranges)
	chrs = [x[0] for x in ranges]
	chrs_types = list(set(chrs))
	for chrs_type in chrs_types: 
		if species == 'hg19':
			f = mmap_fasta("/WPS/BP/yangmingyu/database/chr/"+chrs_type+'.fa')
		if species == 'mm10':
			f = mmap_fasta("/WPS/BP/yangmingyu/database/chr/"+chrs_type+'.fa')
		for id,range in enumerate(ranges):
			if range[0] == chrs_type:
				seqs = f.getseq(range[1]-1,range[2]).upper()
				if range[3] == '-':
					seqs = reverse_complementary(seqs)
				outs[id] = seqs
	return outs

def genome_ranges_2_fa_file(species,a_ranges,file,id_prefix):
	ranges = {}
	ranges_keys = []
	for i,r in enumerate(a_ranges):
		ranges[id_prefix+'#'+str(i)] =  {'chr': r[0], 'right': r[3], 'left': r[2], 'strand': r[1]}
		ranges_keys.append(id_prefix+'#'+str(i))
	get_sequenecs_from_genome(species,ranges)
	f = open(file,'w')
	for k in ranges_keys:
		print >>f,">"+k+"\n"+ranges[k]['seq']
	f.close()

def get_seq_of_range_hg(in1):
	a = re.split(':|-',in1)
	print in1,get_seq_hg_mmap(a[0],int(a[1]),int(a[1])+100)
	print in1,get_seq_hg_mmap(a[0],int(a[2])-100,int(a[2]))

