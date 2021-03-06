import re,copy,sys
import infra00_ranges_operate as in0
import pysam
import d00_sample as d00
import mmap


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

def get_seq_mm_mmap(chr,start,end):
	f = mmap_fasta("/data/Analysis/fanxiaoying/database/mm10/01.bowtie/"+chr+'.fa')
	return f.getseq(start-1,end)

def get_seq_hg_mmap(chr,start,end):
	f = mmap_fasta("/data/Analysis/fanxiaoying/database/hg19/00.genome/"+chr+'.fa')
	return f.getseq(start-1,end).upper()

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

def get_pos_seqs(cursor,specise,ref,poses,server="TANG"):
    refpath = d00.get_ref(cursor,specise,'fa',ref,server)
    ref = pysam.Fastafile(refpath)
    seqs = []
    for pos in poses:
        try:
            seq = ref.fetch(pos[0],pos[1]-1,pos[2]).upper()
        except:
            seq = "NNNN"
            print "Failed To Get Sequence!"
        if pos[3] == "-":
            seq = reverse_complementary(seq)
        seqs.append(seq)
    return seqs

