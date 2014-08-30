import MySQLdb as mb
import re,os

def mm10_ENSG_3UTR(cursor,conn,genename,flank_len):
	cursor.execute("select chr,strand,left_pos,right_pos from mm10_ensembl.UTR3 where genename = %s",([genename]))
	return cursor.fetchall()[0]
