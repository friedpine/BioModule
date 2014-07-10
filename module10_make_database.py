import re
class genes_lists_db(dict):
	print "HEHE"
class exon_database(dict):
	def gtf_exon2bed(self,gtf,bed):
		self['exon_index'] = {}
		self['exon_info'] = {}
		self['transc_info'] = {}
		file = open(gtf)
		bedout = open(bed,'w')
		for line in file:
			array = str.split(line)
			id = array[0]+':'+array[3]+'-'+array[4]
			chr = array[0]
			transc = re.findall("transcript_id \"(.*)\";",line)[0]
			transcript = re.split(r'"',array[9])[1]
			if transc not in self['transc_info']:
				self['transc_info'][transc] = {}
				self['transc_info'][transc]['chr'] = chr
				self['transc_info'][transc]['count'] = 0
				self['transc_info'][transc]['strand'] = array[6] 
				self['transc_info'][transc]['exons_pos'] = []
				if array[6] == '-':
					current = 0
					add_step = -1
				if array[6] == '+':
					current = 1
					add_step = 1
			array[3] = int(array[3])
			array[4] = int(array[4])
			if len(chr)>5:
				continue
			if array[2] == 'exon':
				exon_index = current
				self['transc_info'][transc]['count'] += 1
				self['transc_info'][transc]['exons_pos'].append(array[3])
				self['transc_info'][transc]['exons_pos'].append(array[4])
				if chr not in self['exon_index']:
					self['exon_index'][chr] = {}
					self['exon_info'][chr] = {}
				if array[3] not in self['exon_info'][chr]:
					self['exon_info'][chr][array[3]] = {}
					self['exon_info'][chr][array[3]][array[4]] = {}
					self['exon_info'][chr][array[3]][array[4]]['trans'] = []
					self['exon_info'][chr][array[3]][array[4]]['strand'] = array[6]
					self['exon_info'][chr][array[3]][array[4]]['trans'].append(transcript+'#'+str(exon_index))
					print >>bedout, "%s\t%s\t%s" %(chr,int(array[3])-1,array[4])
					left = int(int(array[3])/1000000)
					if left not in self['exon_index'][chr]:
						self['exon_index'][chr][left] = {}
					self['exon_index'][chr][left][array[3]] = []
					self['exon_index'][chr][left][array[3]].append(array[4])
				elif array[4] not in self['exon_info'][chr][array[3]]:
					self['exon_info'][chr][array[3]][array[4]] = {}
					self['exon_info'][chr][array[3]][array[4]]['trans'] = []
					self['exon_info'][chr][array[3]][array[4]]['strand'] = array[6]
					self['exon_info'][chr][array[3]][array[4]]['trans'].append(transcript+'#'+str(exon_index))
					print >>bedout, "%s\t%s\t%s" %(chr,int(array[3])-1,array[4])
					left = int(int(array[3])/1000000)
					self['exon_index'][chr][left][array[3]].append(array[4])
				else:
					self['exon_info'][chr][array[3]][array[4]]['trans'].append(transcript+'#'+str(exon_index))
				current += add_step
	def gtf_intron2bed(self,gtf,bed,intron_bed_full,length):
		rec_index = 'intron_index'
		rec_info = 'intron_info'
		self[rec_index] = {}
		self[rec_info] = {}
		self['bed_fa_recover_index'] = {}  ##HELP TO MAP READS TO THE RIGHT INTRON
		bedout = open(bed,'w')
		bed_full = open(intron_bed_full,'w')
		for i in self['transc_info']:
			transcript = i
			site = self['transc_info'][i]['exons_pos']
			chr = self['transc_info'][i]['chr']
			strand = self['transc_info'][i]['strand']
			for k in range(1,self['transc_info'][i]['count']):
				index = 2*k-1
				left = site[index]+1
				right = site[index+1]-1
				if right<left or right-left>1000000:
					continue
				if strand == "+":
					exon_left = index
					exon_right = index+1
				elif strand == "-":
					exon_left = self['transc_info'][i]['count']+1-index
					exon_right = self['transc_info'][i]['count']-index
				if chr not in self[rec_index]:
					self[rec_index][chr] = {}
					self[rec_info][chr] = {}
				if left not in self[rec_info][chr]:
					self[rec_info][chr][left] = {}
					self[rec_info][chr][left][right] = {}
					self[rec_info][chr][left][right]['trans'] = []
					self[rec_info][chr][left][right]['strand'] = strand
					self[rec_info][chr][left][right]['trans'].append(transcript+'#'+str(exon_left)+':'+str(exon_right))
					print >>bed_full, "%s\t%s\t%s" %(chr,left-1,right)
					if right-left <= 2*length:
						print >>bedout, "%s\t%s\t%s" %(chr,left-1,right)
						if chr+':'+str(left-1)+'-'+str(right) not in self['bed_fa_recover_index']:
							self['bed_fa_recover_index'][chr+':'+str(left-1)+'-'+str(right)] = []
						self['bed_fa_recover_index'][chr+':'+str(left-1)+'-'+str(right)].append([chr,left,right])
					else:
						print >>bedout, "%s\t%s\t%s" %(chr,left-1,left+length)
						if chr+':'+str(left-1)+'-'+str(left+length) not in self['bed_fa_recover_index']:
							self['bed_fa_recover_index'][chr+':'+str(left-1)+'-'+str(left+length)] = []	
						self['bed_fa_recover_index'][chr+':'+str(left-1)+'-'+str(left+length)].append([chr,left,right])
						print >>bedout, "%s\t%s\t%s" %(chr,right-length,right)
						if chr+':'+str(right-length)+'-'+str(right) not in self['bed_fa_recover_index']:
							self['bed_fa_recover_index'][chr+':'+str(right-length)+'-'+str(right)] = []
						self['bed_fa_recover_index'][chr+':'+str(right-length)+'-'+str(right)].append([chr,left,right])
					left_index = int(left/1000000)
					if left_index not in self[rec_index][chr]:
						self[rec_index][chr][left_index] = {}
						self[rec_index][chr][left_index][left] = []
						self[rec_index][chr][left_index][left].append(right)
				elif right not in self[rec_info][chr][left]:
					self[rec_info][chr][left][right] = {}
					self[rec_info][chr][left][right]['trans'] = []
					self[rec_info][chr][left][right]['strand'] = strand
					self[rec_info][chr][left][right]['trans'].append(transcript+'#'+str(exon_left)+':'+str(exon_right))
					print >>bed_full, "%s\t%s\t%s" %(chr,left-1,right)
					if right-left <= 2*length:
						print >>bedout, "%s\t%s\t%s" %(chr,left-1,right)
						if chr+':'+str(left-1)+'-'+str(right) not in self['bed_fa_recover_index']:
							self['bed_fa_recover_index'][chr+':'+str(left-1)+'-'+str(right)] = []
						self['bed_fa_recover_index'][chr+':'+str(left-1)+'-'+str(right)].append([chr,left,right])
					else:
						print >>bedout, "%s\t%s\t%s" %(chr,left-1,left+length)
						if chr+':'+str(left-1)+'-'+str(left+length) not in self['bed_fa_recover_index']:
							self['bed_fa_recover_index'][chr+':'+str(left-1)+'-'+str(left+length)] = []
						self['bed_fa_recover_index'][chr+':'+str(left-1)+'-'+str(left+length)].append([chr,left,right])
						print >>bedout, "%s\t%s\t%s" %(chr,right-length,right)
						if chr+':'+str(right-length)+'-'+str(right) not in self['bed_fa_recover_index']:
							self['bed_fa_recover_index'][chr+':'+str(right-length)+'-'+str(right)] = []
						self['bed_fa_recover_index'][chr+':'+str(right-length)+'-'+str(right)].append([chr,left,right])
					left_index = int(left/1000000)
					self[rec_index][chr][left_index][left].append(right)
				else:
					self[rec_info][chr][left][right]['trans'].append(transcript+'#'+str(exon_left)+':'+str(exon_right))
	def read_exon_fasta_file(self,fafile):
		file = open(fafile)
		for line in file:
			line = line.strip('\n')
			if re.match('>',line):
				temp = re.split('>|:|-',line)
				chr = temp[1]
				left = int(temp[2])+1
				right = int(temp[3])
			else:
				try:
					self['exon_info'][chr][left][right]['seq'] = line
				except:
					print "NO_recognize",chr,left,right
	def read_intron_fasta_file(self,rec_info,fafile):
		file = open(fafile)
		for line in file:
			line = line.strip('\n')
			if re.match('>',line):
				temp = re.split('>',line)[1]
				continue
			if self['bed_fa_recover_index'][temp] == []:
				print "DELEATED_ALREADY"
			for infos in self['bed_fa_recover_index'][temp]:
				print temp,infos
				if 'seq' not in self[rec_info][infos[0]][infos[1]][infos[2]]:
					self[rec_info][infos[0]][infos[1]][infos[2]]['seq'] = line
				else:
					self[rec_info][infos[0]][infos[1]][infos[2]]['seq'] += ('#'+line)
			self['bed_fa_recover_index'][temp] = []
	def insert_transcript_exon_count_info(self,rec):
		for chr in self[rec]:
			for left in self[rec][chr]:
				for right in self[rec][chr][left]:
					new_transc_info = []
					for t in self[rec][chr][left][right]['trans']:
						transc = re.split('#',t)[0]
						index = int(re.split('#',t)[1])
						count = self['transc_info'][transc]['count']
						new = ''
						if self['transc_info'][transc]['strand'] == '-':
							new = transc+'#'+str(index+count)+'/'+str(count)
						if self['transc_info'][transc]['strand'] == '+':
							new = transc+'#'+str(index)+'/'+str(count)
						new_transc_info.append(new)
					self[rec][chr][left][right]['trans_seq'] = new_transc_info
