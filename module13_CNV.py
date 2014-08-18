import pysam
import os, sys, re
import MySQLdb as mb

conn=mb.connect(host="localhost",user="root",passwd="123456",db='lihang')
cursor = conn.cursor()

chrom_size="hg19_ensembl.hg19"
bamfile="/datc/fanxiaoying/project/project04_LH/01.pipeline/BAMbam_bwa_genome_M10.bam"
tablename = "lihang.CNV_M10_1M"
binsize = 1000000

def GET_BINNED_COUNTS_ALONG_GENOME(cursor,conn,bamfile,binsize,chrom_size,tablename)
  sql_table = "CREATE TABLE "+tablename+" (chr varchar(10),bin_id int, counts int)"	
  cursor.execute("select chr,`length` from "+chrom_size")
  results = cursor.fetchall()
  samfile = pysam.Samfile(bamfile, "rb")
  counts = []
  for chr in results:
    for pos in range(1,int(chr[1]),binsize):
      counts.append((chr[0],int(pos/binsize),samfile.count(chr[0],pos,min(pos+binsize,int(chr[1])))))
  try:
    cursor.execute(sql_table)
  except:
    print "EXISTS"
	cursor.executemany("insert ignore into "+tablename+" values(%s,%s,%s) """,counts)
	conn.commit()
  
      
