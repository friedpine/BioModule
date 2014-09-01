import os,re,sys
import module02_coverage as m02
import d01_geneinfo as d01

def Smooth_By_Windows():
  print "SMOOTH"

def Downhills(cursor,conn,samples,bam_handles,genename,flanksize,min_len,min_sep,rec):
	UTR3 = d01.mm10_refGene_3UTR(cursor,conn,genename,0)
	if UTR3 == ():
		print "NO TRANSC"
		return 0
	chr = UTR3[1][3:]
	print UTR3
	if UTR3[2] == "+":
		print UTR3[3],UTR3[4],UTR3[4]+flanksize
		range_l = UTR3[3]
		range_r = UTR3[4]+flanksize
	elif UTR3[2] == "-":
		range_l = UTR3[3]-flanksize
        range_r = UTR3[4]
	print [chr,range_l,range_r]
	datas = m02.Depth_Data2(samples,bam_handles,[chr,range_l,range_r])
	frames = m02.Depth_Data2_Process_transcript(datas,samples,[range_l,range_r],[],UTR3[2])
	m02.Plot_Depth_Data(samples,frames,rec+genename+'.png')
	m02.Look_For_Downhills(frames,samples,min_len,0,min_sep)
  	
