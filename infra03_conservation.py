import sys
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import infra00_ranges_operate as in0
import subprocess,re
import numpy as np
import infra01_pos2info as in1
import cPickle as pickle
shortNAME = {'Simple_repeat':'SR','Low_complexity':'LC','ERVL-MaLR':'EM','hAT-Charlie':'hC','Y-chromosome':'Yc','TcMar-Tigger':'TT'}

def phastcon_mm10_exon_average(chr,start,end):
	cmd = '/data/Analysis/fanxiaoying/software/bigWigToWig /data/Analysis/fanxiaoying/database/mm10/07.phastcon/mm10.60way.phastCons.bw stdout -chrom=%s -start=%s -end=%s'	%(chr,start,end)
	p1 = subprocess.Popen(cmd,shell='True',stdout=subprocess.PIPE)
	out = []
	for line in p1.stdout:
		if re.match('f',line):
			continue
		out.append(round(float(line),2))
	if out!=[] and round(np.mean(out),2)>0:
		return round(np.mean(out),2)
	else:
		return 0

def repeat_mask_of_genome_ranges(species,ranges):
	events = ranges.keys()
	ranges['index'] = {}
	b_file = 'repeat_mask_task.bed'
	u_file = 'repeat_selected_mask_task.bed'
	f = open(b_file,'w')
	for event in events:
		ranges[event]['repeat'] = {}
		chr = ranges[event]['chr']
		if chr not in ranges['index']:
			ranges['index'][chr] = {}
		left_M_raw = int(ranges[event]['left']/1000000)
		for left_M in [left_M_raw-1,left_M_raw,left_M_raw+1]:
			if left_M not in ranges['index'][chr]:
				ranges['index'][chr][left_M] = []
			ranges['index'][chr][left_M].append(event)
		print >>f,chr+'\t'+str(ranges[event]['left'])+'\t'+str(ranges[event]['right'])+'\n'
	f.close()
	if species == 'mm10':
		bed_file = '/data/Analysis/fanxiaoying/database/mm10/08.repeat/repeat_mask.bed'
	cmd = '/data/Analysis/fanxiaoying/software/bedtools-2.17.0/bin/intersectBed -a %s -b %s -u >%s'	%(bed_file,b_file,u_file)
	subprocess.call(cmd,shell='True')
	file = open(u_file)
	for line in file:
		array = re.split('\s+',line)
		chr = array[0]
		left_M = int(int(array[1])/1000000)
		line_pos = [int(array[1]),int(array[2])]
		for event in ranges['index'][chr][left_M]:
			event_pos = [ranges[event]['left'],ranges[event]['right']]
			info = in0.ranges_overlap(line_pos,event_pos)
			if info[0] == 'T':
				ranges[event]['repeat'][int(array[1])] = array
				#if ranges[event]['repeat'] == {}:
				#	ranges[event]['repeat'][info[2]] = array
				#elif info[2]>ranges[event]['repeat'].keys()[0]:
				#	ranges[event]['repeat'][info[2]] = array
	del ranges['index']
	return ranges

def repeat_mask_of_exon_introns(species,db,segs):
	ranges = {}
	for event in segs:
		poses = in1.get_exons_pos(db,segs[event]['transc'],segs[event]['ids'])
		for id in segs[event]['ids']:
			ranges[event+'#'+str(id)] = {'chr':poses['chr'],'left':poses[id][0],'right':poses[id][1]}
	ranges = repeat_mask_of_genome_ranges(species,ranges)
	for i in ranges:
		if re.split('#',i):
			event = re.split('#',i)[0]
			id = float(re.split('#',i)[1])
			segs[event][id] = {}
			[segs[event][id]['str'],segs[event][id]['raw'],segs[event][id]['ranges'],segs[event][id]['chr']] = simplify_the_repeat_result(ranges[i]['repeat'],shortNAME)
			segs[event][id]['phastcon'] = []
			for r in segs[event][id]['ranges']:
				segs[event][id]['phastcon'].append(phastcon_mm10_exon_average(segs[event][id]['chr'],r[0],r[1]))
	return segs

def simplify_the_repeat_result(repeats,shortNAME):
	simp = []
	ranges = []
	chr = ''
	for left in sorted(repeats.keys()):
		r_name = repeats[left][6]	
		if r_name in shortNAME:
			r_name = shortNAME[r_name]
		chr = repeats[left][0]
		simp.append(r_name+repeats[left][3])
		ranges.append([int(repeats[left][1]),int(repeats[left][2])])
	a = in0.newfunc_set_count_sort(simp)
	return a,simp,ranges,chr

def derive_info_to_ranges_arrays(ranges,ids,all_positive):
	out_all = []
	for id in ids:
		out = []
		for i,name in enumerate(ranges[id]['raw']):
			if all_positive == 'ALL_POS':
				out.append([ranges[id]['chr'],'pos',ranges[id]['ranges'][i][0],ranges[id]['ranges'][i][1]])
			elif all_positive == 'NO':
				if name[-1] == '+':
					out.append([ranges[id]['chr'],'pos',ranges[id]['ranges'][i][0],ranges[id]['ranges'][i][1]])
				else:
					out.append([ranges[id]['chr'],'neg',ranges[id]['ranges'][i][0],ranges[id]['ranges'][i][1]])
		out_all.append(out)
	return out_all

def derive_info_count_repeat_numbers(arrays,repeat_names):
	count = 0
	for a in arrays:
		for n in repeat_names:
			if a[:-1] == n:
				count += 1
	return count
	
def interactions_of_repeats_between_two_sequences(rep1,rep2,considered_repeats):
	if rep1 == '' or rep2 == '' or rep1 == 'NA' or rep2 == 'NA':
		return 0
	elif re.findall('#',rep1) and re.findall('#',rep2):
		out = {}
		a1 = re.split('#',rep1)[:-1]
		a2 = re.split('#',rep2)[:-1]
		a1_d = {}
		a2_d = {}
		for i in a1:
			a1_d[re.split('\*',i)[1]] = int(re.split('\*',i)[0])
		for i in a2:
			a2_d[re.split('\*',i)[1]] = int(re.split('\*',i)[0])
		for i in a1_d:
			if i[-1] == '+':
				mate = i[:-1]+'-'
			if i[-1] == '-':
				mate = i[:-1]+'+'
			if mate in a2_d:
				if i[:-1] not in out:
					out[i[:-1]] = a1_d[i]*a2_d[mate]
				else:
					out[i[:-1]] += a1_d[i]*a2_d[mate]
		out_sort = sorted(out.items(),key=lambda d:d[1],reverse=True)
		t = ''
		count = 0
		for i in out_sort:
			t = t+str(i[1])+'*'+i[0]+'#'
			if i[0] in considered_repeats:
				count += i[1]
		if t == '':
			return 0,0
		else:
			return t,count
	else:
		return 0,0
