import MySQLdb as mb
import re,os

def mm10_refGene_3UTR(cursor,conn,genename,flank_len):
	cursor.execute("select transc,chr,strand,`left`,`right` from mm10_refseq.utr3 where gene = %s",([genename]))
	outs = {}
	for x in cursor.fetchall():
		pos = str(x[3])+'_'+str(x[4])
		if pos in outs:
			continue
		else:
			if x[2] == "+":
				range_flank = [x[3],x[4]+flank_len]
			if x[2] == "-":
				range_flank = [x[3]-flank_len,x[4]]
			outs[pos] = {'chr':x[1],'range':[x[3],x[4]],'range_flank':range_flank,'strand':x[2],'transc':x[0]}
	return outs	
				
def PLOT_APA_Predicted():
	print "HEHE"
