from __future__ import division
import sys
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import pairwise2
from pairwise2 import format_alignment
import subprocess
import cPickle as pickle
import string
import sys,os,re,time

dat_file = sys.argv[1]
exondb = sys.argv[2]
genome_up = int(sys.argv[3])
refseq_up = int(sys.argv[4])
precision = sys.argv[5]
sample = sys.argv[6]

out = pickle.load(open(dat_file))
exondb =  pickle.load(open(exondb))

count_equal = 0
count_unequal = 0
out['event'] = {}
out['event'][sample] = {}
ff = open("t",'a')
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
for chr in out[sample]:
        if chr not in exondb['exon_index']:
                continue
        for left in out[sample][chr]:
		if out[sample][chr][left]['mapp_info'][0]>genome_up or out[sample][chr][left]['mapp_info'][1]>refseq_up:
			continue
                left_M = int(left/1000000)
                right = out[sample][chr][left]['end']
                seql_rev = 'left_NA'
                seqr_rev = 'right_NA'
                exon_l_info = "chrNA"
                exon_r_info = "chrNA"
                left_exon_transc = []
                right_exon_transc = []
                if left_M not in exondb['exon_index'][chr]:
                        continue
                for exl in exondb['exon_index'][chr][left_M]:
                        for exr in exondb['exon_info'][chr][exl]:
                                seq = exondb['exon_info'][chr][exl][exr]['seq']
                                if (left <= exr) and (left > exl) and (left-exl<100):
                                        seql = seq[0:(left-exl+25)]
                                        seql_rev = seql[::-1].upper()
                                        exon_l_info = chr+":"+str(exl)+"-"+str(exr)
                                        LEFT = exl
                                        strand = exondb['exon_info'][chr][exl][exr]['strand']
                                        for tt in exondb['exon_info'][chr][exl][exr]['trans_seq']:
                                                if tt not in left_exon_transc:
                                                        left_exon_transc.append(tt)
                                if (right <exr) and (right >=exl) and (exr-right<100):
                                        seqr = seq[len(seq)-(exr-right+25)-1:len(seq)]
                                        seqr_rev = seqr[::-1].upper()
                                        exon_r_info = chr+":"+str(exl)+"-"+str(exr)
                                        RIGHT = exr
                                        for tt in exondb['exon_info'][chr][exl][exr]['trans_seq']:
                                                if tt not in right_exon_transc:
                                                        right_exon_transc.append(tt)
                seq_observed = out[sample][chr][left]['seq']
                if out[sample][chr][left]['strand'] == 'neg':
                        seq_theory = complementary(seql_rev+seqr_rev)
                if out[sample][chr][left]['strand'] == 'pos':
                        seq_theory = (seql_rev+seqr_rev)[::-1]
                if abs(len(seq_observed)-len(seq_theory))>=6 or (exon_l_info == "chrNA") or (exon_r_info == "chrNA"):
                        continue
		mate_exons = []
		for mate_map in out[sample][chr][left]['mate_genome']:
			a2 = re.split('#',mate_map)
			mate_chr = 'chr'+a2[1]
			if mate_chr != chr:
				print chr,a2,'DIFFERENT_GENOME'
				continue
			pos_left = int(a2[2])
			t01 = re.split('[A-Z]',a2[3])
			t02 = [int(i) for i in t01[0:-1]]
			t03 = [sum(t02[0:i])+pos_left for i in range(len(t02))]
			t04 = re.split('\d+',a2[3])[1:]
			t05 = [i for i in range(len(t04)) if t04[i] =='M']
			for i in t05:
				mate_pos = t03[i]
				mate_pos_M = int(mate_pos/1000000)
				for exon_left in exondb['exon_index'][chr][mate_pos_M]:
					for exon_right in exondb['exon_info'][chr][exon_left]:
						if mate_pos<=exon_right and mate_pos>exon_left-10:
							for tt in exondb['exon_info'][chr][exon_left][exon_right]['trans_seq']:
								if tt not in mate_exons:
									mate_exons.append(tt)
                if precision == 'False':
                        if (exon_l_info != "chrNA") and (exon_r_info != "chrNA") and seq_observed!=seq_theory and abs(len(seq_observed)-len(seq_theory))<6:
                                t0 = time.time()
                                alignment = pairwise2.align.globalms(seq_theory[25:75],seq_observed[25:75],1,-1,-1,-1,one_alignment_only=1)
                                t1 = time.time()
                                score = alignment[0][2]/alignment[0][4]
                                if score <0.8:
                                        continue
                                print "ALIGN_SEQUNECES",t1-t0
                if 1 == 1:
                        print exon_l_info
                        out[sample][chr][left]['event'] = exon_l_info+"_"+exon_r_info
                        out[sample][chr][left]['seq_theory'] = seq_theory
                        if exon_l_info+"_"+exon_r_info not in out['event'][sample]:
                                out['event'][sample][exon_l_info+"_"+exon_r_info] = {}
                                out['event'][sample][exon_l_info+"_"+exon_r_info]['count'] = 0
                                out['event'][sample][exon_l_info+"_"+exon_r_info]['evidences'] = []
                                out['event'][sample][exon_l_info+"_"+exon_r_info]['mate_info_genome'] = []
				out['event'][sample][exon_l_info+"_"+exon_r_info]['mate_info_refseq'] = []
				out['event'][sample][exon_l_info+"_"+exon_r_info]['mate_genome2exon'] = []
			mate_genome = 'NA'
			mate_refseq = 'NA'
                        for temp in out[sample][chr][left]['mate_genome']:
				array = re.split('#',temp)
				if int(array[2])<RIGHT and int(array[2]>LEFT):
					mate_genome = 'YES'
			for temp in out[sample][chr][left]['mate_refseq']:
                                array = re.split('#',temp)
                                for transc in left_exon_transc:
					if re.findall(array[1],transc):
	                                        mate_refseq = 'YES'
			out['event'][sample][exon_l_info+"_"+exon_r_info]['count'] += 1
                        out['event'][sample][exon_l_info+"_"+exon_r_info]['chr'] = chr
                        out['event'][sample][exon_l_info+"_"+exon_r_info]['start'] = LEFT
                        out['event'][sample][exon_l_info+"_"+exon_r_info]['end'] = RIGHT
                        out['event'][sample][exon_l_info+"_"+exon_r_info]['evidences'].append(left)
                        out['event'][sample][exon_l_info+"_"+exon_r_info]['transc'] = [left_exon_transc,right_exon_transc]
			out['event'][sample][exon_l_info+"_"+exon_r_info]['mate_info_genome'].append(mate_genome)
			out['event'][sample][exon_l_info+"_"+exon_r_info]['mate_info_refseq'].append(mate_refseq)
			out['event'][sample][exon_l_info+"_"+exon_r_info]['mate_genome2exon'].append(mate_exons)
for event in out['event'][sample]:
	uni_transc = []
	mate_transc = []
	mate_transc_names = {}
	for i in out['event'][sample][event]['mate_genome2exon']:
		for j in i:
			mate_transc.append(j)
			if not re.split('#|/',j)[0] in mate_transc_names:
				mate_transc_names[re.split('#|/',j)[0]] = []
			mate_transc_names[re.split('#|/',j)[0]].append(re.split('#|/',j)[1])
	print mate_transc_names
	for transc_left in out['event'][sample][event]['transc'][0]:
		for transc_right in out['event'][sample][event]['transc'][1]:
			transc_left_name = re.split('#',transc_left)[0]
			transc_right_name = re.split('#',transc_right)[0]
			if transc_left_name==transc_right_name and transc_left_name in mate_transc_names:
				total_exons = re.split('#|/',transc_left)[2]
				left_index = re.split('#|/',transc_left)[1]
				right_inex = re.split('#|/',transc_right)[1]
				mate = '/'.join(sorted(mate_transc_names[transc_left_name]))	
				transc_uni_ones = transc_left_name+'#'+total_exons+'#'+left_index+'#'+right_inex+'#'+mate
				uni_transc.append(transc_uni_ones)
	out['event'][sample][event]['uni_transc'] = uni_transc
pickle.dump(out,open(sys.argv[7],"wb"),True)
