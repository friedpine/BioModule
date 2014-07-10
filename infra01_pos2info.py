import re,copy,sys
import cPickle as pickle
import subprocess
import infra00_ranges_operate as in0
import mmap
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
ref = pickle.load(open('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules/Ref.all.dat'))

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

def get_infos(db1,db2,poses):
    for chr in poses:
        for pos in poses[chr]:        
            for type in ['exon','intron']:
                poses[chr][pos][type] = {}
                poses[chr][pos][type]['raw'] = []
                poses[chr][pos][type]['transc'] = []
                poses[chr][pos][type]['genes'] = []
                index = type+'_index'
                info = type+'_info'
                pos_M = int(pos/1000000)
		if pos_M not in db2[index][chr]:
			continue
                for left in db2[index][chr][pos_M]:
                    for right in db2[index][chr][pos_M][left]:
                        if pos <= right and pos >= left:
                            temp = copy.deepcopy(db2[info][chr][left][right])
                            temp['Range_left'] = left
                            temp['Range_right'] = right
                            del temp['seq']
                            poses[chr][pos][type]['raw'].append(temp)
                            for trans_info in temp['trans']:
                                transc = re.split('#',trans_info)[0]
                                gene = db1['ref']['t'][transc]
                                if transc not in poses[chr][pos][type]['transc']:
                                    poses[chr][pos][type]['transc'].append(transc)
                                if gene not in poses[chr][pos][type]['genes']:
                                    poses[chr][pos][type]['genes'].append(gene)
                            temp = {}
def get_exons_sequences(db,NM,exon_indexs):
	out = {}
	info = db['transc_info'][NM]
	chr = info['chr']
	strand = info['strand']
	counts = info['count']
	exons_pos = info['exons_pos']
	if strand == '+':
		for index in exon_indexs:
			id = index
			out[index]=db['exon_info'][chr][exons_pos[2*id-2]][exons_pos[2*id-1]]['seq'].upper()
	if strand == '-':
		for index in exon_indexs:
			id = counts+1-index
			out[index]=reverse_complementary(db['exon_info'][chr][exons_pos[2*id-2]][exons_pos[2*id-1]]['seq'].upper())
	return out
def db_get_introns_sequences(refdb,transc,indexes):
	print "hehe"

def get_introns_sequences(db,species,NM,intron_indexs):
	ranges = {}
	info = db['transc_info'][NM]
	chr = info['chr']
        strand = info['strand']
        counts = info['count']
        exons_pos = info['exons_pos']
	out = {}
	if strand == '+':
		for index in intron_indexs:
			id = index
			positions = chr+':'+str(exons_pos[2*id-1]+1)+'-'+str(exons_pos[2*id]-1)
			ranges[positions] = {}
			ranges[positions]['chr'] = chr
			ranges[positions]['left'] = exons_pos[2*id-1]+1
			ranges[positions]['right'] = exons_pos[2*id]-1
			ranges[positions]['strand'] = 'pos'
			ranges[positions]['seq'] = ''
			out[index] = positions
	if strand == '-':
		for index in intron_indexs:
			id = counts-index
			print exons_pos,index,id
			positions = chr+':'+str(exons_pos[2*id-1]+1)+'-'+str(exons_pos[2*id]-1)
			ranges[positions] = {}
			ranges[positions]['chr'] = chr
			ranges[positions]['left'] = exons_pos[2*id-1]+1
			ranges[positions]['right'] = exons_pos[2*id]-1
			ranges[positions]['strand'] = 'neg'
			ranges[positions]['seq'] = ''
			out[index] = positions
	range2 =  get_sequenecs_from_genome(species,ranges)
	for i in out.keys():
		out[i] = range2[out[i]]['seq']
	return out
def get_exons_pos(db,NM,exon_indexs):
	out = {}
	info = db['transc_info'][NM]
        chr = info['chr']
        strand = info['strand']
        counts = info['count']
        exons_pos = info['exons_pos']
	out['chr'] = chr
	for index in exon_indexs:
		if index%1==0:
			if strand == '+':
				id = index
			elif strand == '-':
				id = counts+1-index
			out[index]=[exons_pos[2*id-2],exons_pos[2*id-1]]
		else:
			idx = int(index)
			if strand == '+':
				id = idx
			if strand == '-':
				id = counts-idx
			out[index] = [exons_pos[2*id-1]+1,exons_pos[2*id]-1]
	return out
def get_introns_pos(db,NM,indexs):
	info = db['transc_info'][NM]
	chr = info['chr']
        strand = info['strand']
        counts = info['count']
        exons_pos = info['exons_pos']
        out = {}
	for index_raw in indexs:
		index = int(index_raw)
		if strand == '+':
			id = index
		if strand == '-':
			id = counts-index
		out[index_raw] = [exons_pos[2*id-1]+1,exons_pos[2*id]-1]
	return out
def transform_exon_intron_cordinate_to_genome(db,NM,index,ranges):
	info = db['transc_info'][NM]
	strand = info['strand']
	ranges = in0.ranges_regulate(ranges)
	if index%1==0:
		start_end = get_exons_pos(db,NM,[index])[index]
	else:
		index = int(index)
		start_end = get_introns_pos(db,NM,[index])[index]
	out = []
	if strand == '+':
		out.append(ranges[0]+start_end[0]-1)
		out.append(ranges[1]+start_end[0]-1)
	else:
		out.append(start_end[1]-ranges[1]+1)
		out.append(start_end[1]-ranges[0]+1)
	return out
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

