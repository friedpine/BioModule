import re,copy,sys
import cPickle as pickle
import subprocess
import infra00_ranges_operate as in0
import mmap
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
#ref = pickle.load(open('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/Ref.all.dat'))

class mmap_fasta(object):
    def __init__(self,fname):
        f = file(fname)
        header = f.readline()
        row = f.readline()

        self.ofs = len(header)
        self.lline = len(row)
        self.ldata = len(row.strip())
        self.skip = self.lline-self.ldata
        self.skip_char = row[self.ldata:]
        #print "SKIP",self.skip,self.skip_char
        self.mmap = mmap.mmap(f.fileno(),0,prot=mmap.PROT_READ)

    def getseq(self,start,end):
        l_start = start / self.ldata
        l_end = end / self.ldata
        #print "lines",l_start,l_end
        ofs_start = l_start * self.skip + start + self.ofs
        ofs_end = l_end * self.skip + end + self.ofs
        #print "ofs",ofs_start,ofs_end
        
        s = self.mmap[ofs_start:ofs_end].replace(self.skip_char,"")
        L = end-start
        if len(s) == L:
            return s
        else:
            return s+"N"*(L-len(s))
        return 

def reverse_complementary(seq):
	rc = ''
	reverse = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	for i in seq:
		rc=rc+reverse[i]
	return rc[::-1]

def get_sequenecs_from_genome(species,ranges):
	genome = ref[species]['fa']['genome']
	exon_fa = 'temp_get_sequence.fa'
	bedfile = 'temp_bedfile.bed'
	file = open(bedfile,'w')
	names_pos_id = {}
	for i in ranges.keys():
		print >>file,ranges[i]['chr']+'\t'+str(ranges[i]['left']-1)+'\t'+str(ranges[i]['right'])
		names_pos_id[ranges[i]['chr']+':'+str(ranges[i]['left']-1)+'-'+str(ranges[i]['right'])] = i
	file.close()
	cmd_bed = "/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/fastaFromBed -fi %s -bed %s -fo %s" %(genome,bedfile,exon_fa)
	subprocess.call(cmd_bed,shell=True)
	file = open(exon_fa)
	destination_name = ''
	for line in file:
		if re.match('>',line):
			destination_name = names_pos_id[line[1:-1]]
			ranges[destination_name]['seq'] = ''
		elif ranges[destination_name]['strand'] == 'pos':
			ranges[destination_name]['seq'] = line[:-1].upper()
		elif ranges[destination_name]['strand'] == 'neg':
			ranges[destination_name]['seq'] = reverse_complementary(line[:-1].upper())
	return ranges

def get_sequenecs_from_genome_s(species,ranges):
	genome = ref[species]['fa']['genome']
	exon_fa = '/tmp/temp_get_sequence.fa'
	bedfile = '/tmp/temp_bedfile.bed'
	file = open(bedfile,'w')
	names_pos_id = {}
	for i in ranges:
		print >>file,ranges[i][0]+'\t'+str(ranges[i][2]-1)+'\t'+str(ranges[i][3])
		names_pos_id[ranges[i][0]+':'+str(ranges[i][2]-1)+'-'+str(ranges[i][3])] = i
	file.close()
	cmd_bed = "/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/fastaFromBed -fi %s -bed %s -fo %s" %(genome,bedfile,exon_fa)
	subprocess.call(cmd_bed,shell=True)
	file = open(exon_fa)
	destination_name = ''
	out = {}
	for line in file:
		if re.match('>',line):
			destination_name = names_pos_id[line[1:-1]]
			out[destination_name] = ''
		elif ranges[destination_name][1] == '+':
			out[destination_name] = line[:-1].upper()
		elif ranges[destination_name][1] == '-':
			out[destination_name] = reverse_complementary(line[:-1].upper())
	return out

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
			f = mmap_fasta("/data/Analysis/fanxiaoying/database/hg19/00.genome/"+chrs_type+'.fa')
		if species == 'mm10':
			f = mmap_fasta("/data/Analysis/fanxiaoying/database/mm10/01.bowtie/"+chrs_type+'.fa')
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

