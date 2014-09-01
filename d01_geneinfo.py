import MySQLdb as mb
import re,os

def mm10_refGene_3UTR(cursor,conn,genename,flank_len):
	cursor.execute("select transc,chr,strand,`left`,`right` from mm10_refseq.utr3 where gene = %s",([genename]))
	return cursor.fetchall()[0]
		
