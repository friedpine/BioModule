import sys,os
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import random as rd
import numpy as np
#import module03_rna_procerss as m03
import infra00_ranges_operate as in0
import infra01_pos2info as in1
import infra02_blast_info_proc as in2
import infra03_conservation as in3
import MySQLdb as mb

def random_build_fake_cRNA(refdb,dbname,table1,numbers):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=dbname)
	cursor = conn.cursor()
	cursor.execute("SELECT id,exon_count FROM "+refdb+".transc WHERE exon_count>5 AND id NOT IN (SELECT transc FROM mm_cRNA.01_all_events)")
	r0 = list(cursor.fetchall())
	transc_3k = rd.sample(r0,numbers)
	values = []
	for transc in transc_3k:
		r1 = rd.sample(range(2,transc[1]),2)
		values.append([transc[0],min(r1),max(r1)])
	cursor.executemany("insert into "+table1+" (id,transc,ex_left_id,ex_right_id)values(NULL,%s,%s,%s)",values)
	conn.commit()
def blast_introns_sequence(refdb,table_in,dbname,table1,table2,evalue,wordsize):
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=dbname)
	cursor = conn.cursor()
	cursor.execute("select * from "+table1)
	r0 = cursor.fetchall()
	for t in r0:
		print '###########',t
		try:
			ranges = {}
			cursor.execute("select * from "+refdb+'.'+table_in+" where transc = %s and ex_order = %s ",[t[1],t[2]-1])
			r1 = cursor.fetchall()[0]
			ranges['r2'] = [r1[2],'+',r1[4],r1[5]]
			cursor.execute("select * from "+refdb+'.'+table_in+" where transc = %s and ex_order = %s ",[t[1],t[3]])
			r1 = cursor.fetchall()[0]
			ranges['r1'] = [r1[2],'+',r1[4],r1[5]]
			in2.blast_genome_positions('mm10',ranges,evalue,wordsize,'RC')
			result = in2.blast_fmt7_out_read_for_db('./temp_blast_r1r2.txt',40,'RC')
			if result != []:
				for r2 in result:
					cursor.execute("insert into "+table2+" values(%s,%s,%s,%s,%s,%s,%s,%s,%s)",[t[0]]+r2)
		except:
			cursor.execute("insert into "+table2+" values(%s,0,0,0,0,0,0,0,0)",[t[0]])
		conn.commit()
# class circ_explore(m03.cirrna_meth1):
# 	def __init__(self,input):
# 		for key_used in input.keys():
# 			print key_used
# 			self[key_used] = input[key_used]
# 		print "Welcome to m03_sub_class"
# 	def blast_exons_introns_of_transcript(self,events_rec,exon_db,species,folder,evalue,wordsize,file_rec,rec):
# 		if not os.path.exists(folder):
# 			os.mkdir(folder)
# 		self[rec] = []
# 		for event in self[events_rec]:
# 			fasta = folder+'/'+event+'.fa'
# 			db_file = folder+'/'+event+'.db'
# 			out = folder+'/'+file_rec+'_'+event+'.txt'
# 			self[events_rec][event][file_rec] = out
# 			if not os.path.exists(fasta):
# 				p_transc = self[events_rec][event]['transc']['prefered']
# 				total_exons = self[events_rec][event]['transc'][p_transc]['total']
# 				start = self[events_rec][event]['transc'][p_transc]['start']
# 				end = self[events_rec][event]['transc'][p_transc]['end']
# 				print p_transc,range(1,total_exons+1)
# 				sequences = in1.get_exons_sequences(exon_db,p_transc,range(max(1,start-1),min(total_exons+1,end+2)))
# 				f = open(fasta,'w')
# 				for exon_id in range(max(1,start-1),min(total_exons+1,end+2)):
# 					f.write('>'+str(total_exons)+"#"+str(exon_id)+'\n'+sequences[exon_id]+'\n')
# 				sequences = in1.get_introns_sequences(exon_db,species,p_transc,range(max(1,start-1),min(total_exons,end+1)))
# 				for intron_id in range(max(1,start-1),min(total_exons,end+1)):
# 					f.write('>'+str(total_exons)+"#"+str(intron_id+0.5)+'\n'+sequences[intron_id]+'\n') 
# 				f.close()
# 			cmd1 = "/data/Analysis/fanxiaoying/software/ncbi-blast-2.2.28+/bin/makeblastdb -dbtype nucl -in %s -out %s " %(fasta,db_file)
# 			if not os.path.exists(db_file+'.nhr'):
# 				self[rec].append(cmd1)
# 			cmd2 = "/data/Analysis/fanxiaoying/software/ncbi-blast-2.2.28+/bin/blastn -query %s -task blastn -db %s -outfmt 7 -evalue %s -word_size %s -out %s " %(fasta,db_file,evalue,wordsize,out)
# 			if not os.path.exists(out):
# 				self[rec].append(cmd2)
# 	def get_A012_B012_repeat_mask(self,db,rec):
# 		segs = {}
# 		for event in self['event_info']:
# 			self['event_info'][event][rec] = {}
# 			transc = self['event_info'][event]['transc']['prefered']
# 			info = self['event_info'][event]['transc'][transc]
# 			if info['start']==1 or info['end']==info['total']:
# 				continue
# 			indexs = [info['start']-0.5,info['start'],info['start']+0.5,info['end']-0.5,info['end'],info['end']+0.5]
# 			segs[event] = {'transc': transc, 'ids':indexs}
# 		ranges = in3.repeat_mask_of_exon_introns('mm10',db,segs)
# 		for event in ranges:
# 			self['event_info'][event][rec] = ranges[event]
# 	def test_nine_models(self,rec,file_rec,samples_limits,reads_limits):
# 		self[rec] = [[[],[],[],[],[]],[[],[],[],[],[]],[[],[],[],[],[]],[[],[],[],[],[]],[[],[],[],[],[]]]
# 		self[rec+'_count'] = [[[],[],[],[],[]],[[],[],[],[],[]],[[],[],[],[],[]],[[],[],[],[],[]],[[],[],[],[],[]]]
# 		count = 0
# 		for event in self['event_info']:
# 			self['event_info'][event][rec] = {}
# 			if len([i for i in self['event_info'][event]['count'] if i>=reads_limits]) < samples_limits:
# 				continue
# 			info = self['event_info'][event]['transc'][self['event_info'][event]['transc']['prefered']]
# 			As_exon = [max(info['start']-1,1),info['start'],min(info['start']+1,info['total'])]
# 			Bs_exon = [max(info['end']-1,1),info['end'],min(info['end']+1,info['total'])]
# 			if Bs_exon[0] <= As_exon[2] or len(set(As_exon))==2 or len(set(Bs_exon))==2:
# 				continue
# 			count += 1
# 			As = [As_exon[0],As_exon[0]+0.5,As_exon[1],As_exon[1]+0.5,As_exon[2]]
# 			Bs = [Bs_exon[0],Bs_exon[0]+0.5,Bs_exon[1],Bs_exon[1]+0.5,Bs_exon[2]]
# 			print event,As,Bs
# 			all = in2.exons(self['event_info'][event][file_rec],info['total'])
# 			for i in range(len(As)):
# 				for j in range(len(Bs)):
# 					value = [0]
# 					A_range = []
# 					B_range = []
# 					Ai = As[i]
# 					Bj = Bs[j]
# 					seq = [0]*all[Ai]['length']
# 					if Ai in all:
# 						for index,type in enumerate(all[Ai]['type']):
# 							if type == 'RC' and all[Ai]['target'][index] == Bj:
# 								left = min(all[Ai]['qrange'][index][0],all[Ai]['qrange'][index][1])
# 								right = max(all[Ai]['qrange'][index][0],all[Ai]['qrange'][index][1])
# 								left_intron = min(all[Ai]['srange'][index][0],all[Ai]['srange'][index][1])
# 								right_intron = max(all[Ai]['srange'][index][0],all[Ai]['srange'][index][1])
# 								seq[left:right] = [1]*(right-left)
# 								value.append(right-left)
# 								B_range.append([left_intron,right_intron])
# 								A_range.append([left,right])	
# 					if Ai != Bj:
# 						self[rec][i][j].append(float(sum(seq))/all[Ai]['length'])
# 						self[rec+'_count'][i][j].append(max(value))
# 					if i == 1 and j == 3:
# 						self['event_info'][event][rec]['model_A2B4'] = {}
# 						self['event_info'][event][rec]['model_A2B4']['A_range'] = A_range
# 						self['event_info'][event][rec]['model_A2B4']['B_range'] = B_range
# 					if i == 1 and j == 1:
# 						self['event_info'][event][rec]['model_A2B2'] = {}
# 						self['event_info'][event][rec]['model_A2B2']['A_range'] = A_range
# 						self['event_info'][event][rec]['model_A2B2']['B_range'] = B_range
# 					if i == 3 and j == 3:
# 						self['event_info'][event][rec]['model_A4B4'] = {}
# 						self['event_info'][event][rec]['model_A4B4']['A_range'] = A_range
# 						self['event_info'][event][rec]['model_A4B4']['B_range'] = B_range
# 					if i == 3 and j == 1:
# 						self['event_info'][event][rec]['model_A4B2'] = {}
# 						self['event_info'][event][rec]['model_A4B2']['A_range'] = A_range
# 						self['event_info'][event][rec]['model_A4B2']['B_range'] = B_range
# 		print "THE_TOTAL_EVENT",count
# 	def model_A0b0_by_single_exon_cRNA(self,rec,file_rec,samples_limits,reads_limits):
# 		self[rec] = [[[0],[0],[0]],[[0],[0],[0]],[[0],[0],[0]]]
# 		self[rec+'_count'] = [[[0],[0],[0]],[[0],[0],[0]],[[0],[0],[0]]]
# 		for event in self['event_info']:
# 			if len([i for i in self['event_info'][event]['count'] if i>=reads_limits]) < samples_limits:
# 				continue
# 			info = self['event_info'][event]['transc'][self['event_info'][event]['transc']['prefered']]
# 			if info['start'] != info['end'] or info['start'] == 1 or info['start'] == info['total']:
# 				continue
# 			all = in2.exons(self['event_info'][event][file_rec],info['total'])
# 			As = [info['start']-0.5,info['start'],info['start']+0.5]
# 			Bs = [info['start']-0.5,info['start'],info['start']+0.5]
# 			print event,As,Bs
# 			all = in2.exons(self['event_info'][event][file_rec],info['total'])
# 			for i in range(len(As)):
# 				for j in range(len(Bs)):
# 					value = [0]
# 					Ai = As[i]
# 					Bj = Bs[j]
# 					seq = [0]*all[Ai]['length']
# 					print self['event_info'][event][file_rec],all[Ai]['length'],Ai,Bj
# 					if Ai in all:
# 						for index,type in enumerate(all[Ai]['type']):
# 							if type == 'RC' and all[Ai]['target'][index] == Bj:
# 								left = min(all[Ai]['qrange'][index][0],all[Ai]['qrange'][index][1])
# 								right = max(all[Ai]['qrange'][index][0],all[Ai]['qrange'][index][1])
# 								seq[left:right] = [1]*(right-left)
# 								value.append(right-left)
# 					if Ai != Bj:
# 						self[rec][i][j].append(float(sum(seq))/all[Ai]['length'])
# 						self[rec+'_count'][i][j].append(max(value))
# 	def model_A2B4_conservation_mask(self,db,model,rec,result):
# 		self[result] = {}
# 		self[result]['signal_A'] = []
# 		self[result]['signal_B'] = []
# 		self['ranges_for_mask'] = {}
# 		for event in self['event_info']:
# 			if model not in self['event_info'][event][rec] or self['event_info'][event][rec][model]['A_range'] == []:
# 				continue
# 			transc = self['event_info'][event]['transc']['prefered']
# 			info = self['event_info'][event]['transc'][transc]
# 			chr = self['event_info'][event]['info']['chr']
# 			[A_sort,id] = in0.get_longest_range(self['event_info'][event][rec][model]['A_range'])
# 			B_sort	= self['event_info'][event][rec][model]['B_range'][id]
# 			A_test_genome = in1.transform_exon_intron_cordinate_to_genome(db,transc,info['start']-0.5,A_sort)	
# 			B_test_genome = in1.transform_exon_intron_cordinate_to_genome(db,transc,info['end']+0.5,B_sort)
# 			self['ranges_for_mask'][event+'#A'] = {'chr': chr, 'right': A_test_genome[1], 'left': A_test_genome[0]}
# 			self['ranges_for_mask'][event+'#B'] = {'chr': chr, 'right': B_test_genome[1], 'left': B_test_genome[0]}
# 			#blast_AB = in2.blast_genome_positions('mm10',[chr,chr],A_test_genome,B_test_genome,10,10,'RC') 
# 			print self['event_info'][event]['exint_e10w10']
# 			a = in3.phastcon_mm10_exon_average(chr,A_test_genome[0],A_test_genome[1])
# 			if a>=0:
# 				self[result]['signal_A'].append(a)
# 			b = in3.phastcon_mm10_exon_average(chr,B_test_genome[0],B_test_genome[1])
# 			if b>=0:
# 				self[result]['signal_B'].append(b)
# 		#print self['ranges_for_mask']
# 		#self['A2B4_repeat'] = in3.repeat_mask_of_genome_ranges('mm10',self['ranges_for_mask'])
# 		f = open('repeat_A2B4.txt','w')
# 		for event in self['event_info']:
# 			if model not in self['event_info'][event][rec] or self['event_info'][event][rec][model]['A_range'] == []:
# 				continue
# 			print >>f,event,
# 			if self['A2B4_repeat'][event+'#A']['repeat'] != {}:
# 				one_key = self['A2B4_repeat'][event+'#A']['repeat'].keys()[0]
# 				info = self['A2B4_repeat'][event+'#A']
# 				info2 = self['A2B4_repeat'][event+'#A']['repeat'][one_key]
# 				print >>f,info['left'],info['right'],info2[1],info2[2],info2[6],in0.ranges_overlap([info['left'],info['right']],[info2[1],info2[2]])[1],
# 				A_test1 = [int(info2[1]),int(info2[2])]
# 				A_test2 = [info['left'],info['right']]
# 			else:
# 				print >>f,'NA',
# 			if self['A2B4_repeat'][event+'#B']['repeat'] != {}: 
# 				one_key = self['A2B4_repeat'][event+'#B']['repeat'].keys()[0]
# 				info = self['A2B4_repeat'][event+'#B']
# 				info2 = self['A2B4_repeat'][event+'#B']['repeat'][one_key]
# 				print >>f,info['left'],info['right'],info2[1],info2[2],info2[6],in0.ranges_overlap([info['left'],info['right']],[info2[1],info2[2]])[1],
# 				B_test1 = [int(info2[1]),int(info2[2])]
# 				B_test2 = [info['left'],info['right']]
# 			else:
# 				print >>f,'NA'
# 			if self['A2B4_repeat'][event+'#A']['repeat'] != {} and self['A2B4_repeat'][event+'#B']['repeat'] != {}:
# 				print 'mm10',[info2[0],info2[0]],['pos','pos'],A_test1,B_test1,A_test2,B_test2,10,10,'RC'
# 				try:
# 					blast_AB_1 = in2.blast_genome_positions('mm10',[info2[0],info2[0]],['pos','pos'],A_test1,B_test1,10,10,'RC')
# 					blast_AB_2 = in2.blast_genome_positions('mm10',[info2[0],info2[0]],['pos','pos'],A_test2,B_test2,10,10,'RC')
# 					print >>f,blast_AB_2[0][1],blast_AB_1[0][1]
# 				except:
# 					print >>f,'NA','NA'
# 			else:
# 				print >>f,'NA','NA'
# 	def repeat_similarity_research(self,all_keys,reads_limits,count_limits,inrec,rec1,rec2,name_space):
# 		for event in all_keys:
# 			transc = self['event_info'][event]['transc']['prefered']
# 			info = self['event_info'][event]['transc'][transc]
# 			if len([i for i in self['event_info'][event]['count'] if i>=reads_limits]) < count_limits:
# 				continue
# 			if info['start']==1 or info['end']==info['total']:
# 				continue
# 			print event,info['start'],info['end'],
# 			self['event_info'][event][rec1] = []
# 			self['event_info'][event][rec2] = []
# 			in_1 = self['event_info'][event][inrec][info['start']-0.5]
# 			in_2 = self['event_info'][event][inrec][info['end']+0.5]
# 			for idx1,r1 in enumerate(in_1['raw']):
# 				for idx2,r2 in enumerate(in_2['raw']):
# 					if r2[:-1] == r1[:-1] and r2[-1] != r1[-1] and r1[:-1] in name_space:
# 						out = [r2[:-1],idx1,idx2,in_1['phastcon'][idx1],in_2['phastcon'][idx2],0]
# 						self['event_info'][event][rec1].append(out)
# 			poses =  in3.derive_info_to_ranges_arrays(self['event_info'][event][inrec],[info['start']-0.5,info['end']+0.5],'ALL_POS')
# 			result =  in2.blast_genome_multi_positions('mm10',poses[1],poses[0],10,10,'RC')
# 			for index,cond in enumerate(self['event_info'][event][rec1]):
# 				index1 = 'r2#'+str(cond[1])
# 				index2 = 'r1#'+str(cond[2])
# 				if index1 not in result:
# 					continue
# 				for target_id,target in enumerate(result[index1]['target']):
# 					if target == index2 and result[index1]['match_len'][target_id]>cond[-1]:
# 						cond[-1] = result[index1]['match_len'][target_id]
# 				self['event_info'][event][rec1][index] = cond
# 			in_1 = self['event_info'][event][inrec][info['start']+0.5]
# 			in_2 = self['event_info'][event][inrec][info['end']-0.5]
# 			for idx1,r1 in enumerate(in_1['raw']):
# 				for idx2,r2 in enumerate(in_2['raw']):
# 					if r2[:-1] == r1[:-1] and r2[-1] != r1[-1] and r1[:-1] in name_space:
# 						out = [r2[:-1],idx1,idx2,in_1['phastcon'][idx1],in_2['phastcon'][idx2],0]
# 						self['event_info'][event][rec2].append(out)
# 			poses =  in3.derive_info_to_ranges_arrays(self['event_info'][event][inrec],[info['start']+0.5,info['end']-0.5],'ALL_POS')
# 			result =  in2.blast_genome_multi_positions('mm10',poses[1],poses[0],10,10,'RC')
# 			for index,cond in enumerate(self['event_info'][event][rec2]):
# 				index1 = 'r2#'+str(cond[1])
# 				index2 = 'r1#'+str(cond[2])
# 				if index1 not in result:
# 					continue
# 				for target_id,target in enumerate(result[index1]['target']):
# 					if target == index2 and result[index1]['match_len'][target_id]>cond[-1]:
# 						cond[-1] = result[index1]['match_len'][target_id]
# 				self['event_info'][event][rec2][index] = cond
# 			print in3.interactions_of_repeats_between_two_sequences(self['event_info'][event][inrec][info['start']-0.5]['str'],self['event_info'][event][inrec][info['end']+0.5]['str'],name_space)
# 			print self['event_info'][event][rec1]
# 			print in3.interactions_of_repeats_between_two_sequences(self['event_info'][event][inrec][info['start']+0.5]['str'],self['event_info'][event][inrec][info['end']-0.5]['str'],name_space)
# 			print self['event_info'][event][rec2]
			
						


