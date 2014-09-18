import os,re,sys
import module02_coverage as m02
import d01_geneinfo as d01
import d02_db_tables as d02
import infra00_ranges_operate as in00
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import copy
import infra01_pos2info as in1

def Smooth_By_Windows():
  print "SMOOTH"

def Look_For_Downhills(datas,samples,min_len,min_depth,mergeable_gap,min_ratio):
	outs = {}
	for sample in samples:
		outs[sample] = []
		d = datas[sample]
		strs  = "".join([str((sum(d[x-10:x])>=sum(d[x:x+10]) and sum(d[x-10:x])>=10)*1) for x in range(10,len(d)-10)])
		strs = '0000000000'+strs
		down_hills = []
		for m in re.finditer('1{1,}',strs):
			if m.end()-m.start()>min_len and sum(d[m.start()-10:m.start()])*min_ratio>sum(d[m.end()-10:m.end()]):
				down_hills.append([m.start(),m.end()])
		#MERGE_THE_RANGES:GAP_LENGTH<mergeable_gap and with near coverage depth
		if down_hills == []:
			continue
		down_hills_merged = []
		merge_infos = '0'
		for x in range(1,len(down_hills)):
			temp = '0'
			down_range = down_hills[x]
			down_range_front = down_hills[x-1]
			if down_range[0]<down_range_front[1]+mergeable_gap:
				depth_range = sum(d[down_range[0]:down_range[0]+10])
				depth_range_front = sum(d[down_range_front[1]-10:down_range_front[1]])
				if abs(depth_range-depth_range_front)<10 or min(depth_range,depth_range_front)> 0.75*max(depth_range,depth_range_front): 
					temp = '1'
			merge_infos = merge_infos+temp
		for m in re.finditer('01{0,}',merge_infos):
			down_hills_merged.append([down_hills[m.start()][0],down_hills[m.end()-1][1]])	
		for i in down_hills_merged:
			outs[sample].append(i)
	return outs

def Fit_For_APA_Site(samples,datas,downhills,min_fit_len):
	APAs = {}
	for sample in samples:
		depth = datas[sample]
		APA = {}
		for down in downhills[sample]:
			if depth[down[1]] == 0:
				down[1] = down[1]-10
				APA[down[1]] = {'pos':down,'depth':[depth[down[0]],0],'color':['cyan','red']}
				APA[down[1]]['area'] = sum(depth[down[0]:down[1]])/(down[1]-down[0])
			else:
				fit_range_left = int(down[0]+(down[1]-down[0])*0.1)
				fit_range_right = int(max(fit_range_left+min_fit_len,down[0]+(down[1]-down[0])*0.9))
				fit_data_x = range(fit_range_left,fit_range_right)
				fit_data_y = depth[fit_range_left:fit_range_right]
				fit_result = np.polyfit(fit_data_x,fit_data_y,1)
				fited_APA_site = int(fit_result[1]/abs(fit_result[0]))
				down[1] = fit_range_right
				down[0] = fit_range_left 
				APA[fited_APA_site] = {'pos':down+[fited_APA_site],'depth':[depth[down[0]],depth[down[1]]]+[0],'color':['cyan','cyan','red']}
				APA[fited_APA_site]['area'] = sum(depth[down[0]:down[1]])/(down[1]-down[0])
		APAs[sample] = APA
	return APAs

def Plot_APA_Downhills(samples,datas,downhills,APAs,points,ymax,filename):
	plt.figure(figsize=(8, 6))
	n = len(samples)
	scale_size = int(len(datas[samples[0]])/points)+1
	ids = [i for i in range(0,len(datas[samples[0]]),scale_size)]
	poses = [datas['pos'][i] for i in ids]
	for number,sample in enumerate(samples):
		depths = datas[sample]
		APA = APAs[sample]
		ax = plt.subplot(n,1,number+1)
		blue_x = []
		blue_y = []
		for downhill in downhills[sample]:
			for site in range(downhill[0],downhill[1]):
				if site in ids:
					blue_x.append(ids.index(site))
					blue_y.append(depths[site])
		gray_x = [x for x in range(0,len(ids)) if x not in blue_x]
		gray_y = [depths[ids[x]]  for x in gray_x]
		ax.bar(gray_x,gray_y,color='gray',edgecolor='gray')
		ax.bar(blue_x,blue_y,color='blue',edgecolor='blue')
		if ymax=="MAX":
			ymax = max(gray_y)
		APA_x = []
		APA_y = []
		APA_color = []
		for site in APA:
			APA_x = APA_x+APA[site]['pos']
			APA_y = APA_y+APA[site]['depth']
			APA_color = APA_color+APA[site]['color']
		APA_x = np.array(APA_x)/scale_size
		#ax.scatter(APA_x,APA_y,c=APA_color,marker=u'*')
		ax.bar(APA_x,APA_y,color=APA_color,edgecolor=APA_color,width=1)
		ax.set_xlim([0,len(ids)])
		ax.set_ylim([0,ymax])
	plt.savefig(filename)
	plt.clf()

def Create_APA_Table(cursor,conn,tablename):
	sql = "CREATE TABLE "+tablename+""" (
  `gene` varchar(30) DEFAULT NULL,
  `transc` varchar(20) DEFAULT NULL,
  `strand` varchar(2) DEFAULT NULL,
  `sample` varchar(20) DEFAULT NULL,
  `chr` varchar(20) DEFAULT NULL,
  `pos` int(11) DEFAULT NULL,
  `down_left` int(11) DEFAULT NULL,
  `down_right` int(11) DEFAULT NULL,
  `down_left_depth` int(11) DEFAULT NULL,
  `down_right_depth` int(11) DEFAULT NULL,
  `mean` int(11) DEFAULT NULL,
  UNIQUE KEY `u` (`sample`,`chr`,`pos`)) ENGINE=InnoDB DEFAULT CHARSET=latin1"""
	try:
		cursor.execute(sql)	
	except:
		print tablename,"EXISTS"

def Save_APAs_info_to_Database(cursor,conn,tablename,genename,frames,utr,samples,APAs):
	for sample in samples:
		for site in APAs[sample]:
			if site>len(frames['pos']):
				continue
			APA = APAs[sample][site]
			APA_pos = frames['pos'][site]
			down_left = frames['pos'][APA['pos'][0]]
			down_right = frames['pos'][APA['pos'][1]]
			cursor.execute("insert ignore into "+tablename+" values(%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)",[genename,utr['transc'],utr['strand'],sample,utr['chr'],APA_pos,down_left,down_right,APA['depth'][0],APA['depth'][1],APA['area']])
	conn.commit()
	
def APAs_Site_Clustering(cursor,conn,sourcetable,outtable,window_size,min_depth,min_supp):
	try:
		cursor.execute("create table "+outtable+"""
			 (	`id` int(11) NOT NULL AUTO_INCREMENT,
			 	`gene` varchar(50) DEFAULT NULL,
				`chr` varchar(20) DEFAULT NULL,
				`strand` varchar(10) DEFAULT NULL,
				`pos` int(11) DEFAULT NULL,
				`sample_count` int(11) DEFAULT NULL,
				`pos_std` int(11) DEFAULT NULL
				 UNIQUE KEY `gene` (`gene`,`pos`),
				KEY `id` (`id`)
				) ENGINE=InnoDB DEFAULT CHARSET=latin1""")
	except:
		print "EXISTS"
	cursor.execute("select * from "+sourcetable)
	All_sites = cursor.fetchall()
	Genes_APAs = {}
	Genes_infos = {}
	#Make a small index for each gene
	for i,info in enumerate(All_sites):
		gene = info[0]
		if info[10]>=min_depth:
			if gene not in Genes_APAs:
				Genes_APAs[gene] = []
				Genes_infos[gene] = [info[4],info[2]]
			Genes_APAs[gene].append(i)
	outs = []
	for gene in Genes_APAs:
		ids = Genes_APAs[gene]
		if len(ids) < min_supp:
			continue
		poses = [All_sites[x][5] for x in ids]
		pos_clusters = in00.clustering_by_windowSize(poses,window_size)
		#print  pos_clusters
		for cluster in pos_clusters:
			if len(cluster)>=min_supp:
				out = [gene]+Genes_infos[gene]+[int(np.median(cluster)),len(cluster),int(np.std(cluster))]
				outs.append(out)
	cursor.executemany("insert into "+outtable+" values(%s,%s,%s,%s,%s,%s)",outs)
	conn.commit()
	
def APAs_Sites_Flanking_Sequences(cursor,conn,tablename,species,flanksize,colname,colinfo):
	cursor.execute("select id,chr,pos,strand from "+tablename)
	APA_sites = cursor.fetchall()
	poses = []
	ids = []
	for site in APA_sites:
		poses.append([site[1],site[2]-flanksize,site[2]+flanksize,site[3]])
		ids.append(site[0])
	seqs = in1.get_multiseqs_mmap(species,poses)	
	d02.append_colume_info_to_tables(cursor,conn,tablename,colname,colinfo,ids,seqs)

def Downhills(cursor,conn,tablename,samples,bam_handles,genenames,min_sample_size,flanksize,min_len,merge_sep,min_ratio,min_fit_len,points,ymax,rec,chr_phase):
	for genename in genenames:
		UTR3 = d01.hg19_refGene_3UTR(cursor,conn,genename,flanksize)
		print genename,UTR3
		if UTR3 == {}:
			print "NO TRANSC"
			return 0
		for pos in UTR3:
			utr = UTR3[pos]
			try:
				read_counts = m02.Depth_Read_Counts(samples,bam_handles,[utr['chr'][chr_phase:]]+utr['range_flank'])
			except:
				continue
			if len([x for x in read_counts if x>=10])<min_sample_size:
				continue
			datas = m02.Depth_Data2(samples,bam_handles,[utr['chr'][chr_phase:]]+utr['range_flank'])
			frames = m02.Depth_Data2_Process_transcript(datas,samples,utr['range_flank'],[],utr['strand'])
			downhills = Look_For_Downhills(frames,samples,min_len,0,merge_sep,min_ratio)
			downhills_c =copy.deepcopy(downhills)
			APAs = Fit_For_APA_Site(samples,frames,downhills,min_fit_len)
			Save_APAs_info_to_Database(cursor,conn,tablename,genename,frames,utr,samples,APAs)
		#Plot_APA_Downhills(samples,frames,downhills_c,APAs,points,ymax,rec+genename+'_'+utr['transc']+'.png')
 
def Downhills_Intermediate(cursor,conn,samples,bam_handles,genename,flanksize,min_len,merge_sep,min_ratio):
	UTR3 = d01.mm10_refGene_3UTR(cursor,conn,genename,flanksize)
	if UTR3 == {}:
		print "NO TRANSC"
		return 0
	for pos in UTR3:
		utr = UTR3[pos]
		datas = m02.Depth_Data2(samples,bam_handles,[utr['chr'][3:]]+utr['range_flank'])
		frames = m02.Depth_Data2_Process_transcript(datas,samples,utr['range_flank'],[],utr['strand'])
		downhills = Look_For_Downhills(frames,samples,min_len,0,merge_sep,min_ratio)
		print downhills
