import os,re,sys
import module02_coverage as m02
import d01_geneinfo as d01
import matplotlib
import matplotlib.pyplot as plt


def Smooth_By_Windows():
  print "SMOOTH"

def Look_For_Downhills(datas,samples,min_len,min_depth,min_seq):
	outs = {}
	for sample in samples:
		outs[sample] = []
        d = datas[sample]
        strs  = "".join([str((sum(d[x-10:x])>=sum(d[x:x+10]) and sum(d[x-10:x])>=10)*1) for x in range(10,len(d)-10)])
        strs = '0000000000'+strs
        down_hills = []
        for m in re.finditer('1{1,}',strs):
            if m.end()-m.start()>min_len:
                down_hills.append([m.start(),m.end()])
		down_hills_sep = [down_hills[i] for i in range(0,len(down_hills)-1) if down_hills[i+1][0]>=down_hills[i][1]+min_seq]    
		down_hills_sep.append(down_hills[-1])
		print down_hills_sep
        for i in down_hills_sep:
			outs[sample].append(i)
			print sample,i
	print outs
	return outs

def Plot_APA_Downhills(samples,datas,downhills,points,ymax,filename):
	plt.figure(figsize=(5, 4))
	n = len(samples)
	scale_size = int(len(datas[samples[0]])/points)+1
	ids = [i for i in range(0,len(datas[samples[0]]),scale_size)]
	poses = [samples['pos'][i] for i in ids]
	for number,sample in enumerate(samples):
		ax = plt.subplot(n,1,number+1)
		plot_y = [datas[sample][x] for x in ids]
		bars = ax.bar(range(0,len(plot_y)),plot_y,label=sample,color='gray')
		for downhill in downhills[sample]:
			for site in range(downhill[0],downhill[1]):
				if site in poses:
					bars[poses.index(site)].set_facecolor('blue')
		ax.set_ylim([0,ymax])
	plt.savefig(filename)
	plt.clf()	


def Downhills(cursor,conn,samples,bam_handles,genename,flanksize,min_len,min_sep,points,ymax,rec):
	UTR3 = d01.mm10_refGene_3UTR(cursor,conn,genename,flanksize)
	if UTR3 == {}:
		print "NO TRANSC"
		return 0
	for pos in UTR3:
		utr = UTR3[pos]
		datas = m02.Depth_Data2(samples,bam_handles,[utr['chr'][3:]]+utr['range_flank'])
		frames = m02.Depth_Data2_Process_transcript(datas,samples,utr['range_flank'],[],utr['strand'])
		downhills = Look_For_Downhills(frames,samples,min_len,0,min_sep)
		#m02.Plot_Depth_Data(samples,frames,rec+genename+'_'+utr['transc']+'.png')
		Plot_APA_Downhills(samples,frames,downhills,points,ymax,rec+genename+'_'+utr['transc']+'.png')
  
	
