import os, sys, re
import MySQLdb as mb
sys.path.insert(0,'/datc/fanxiaoying/git/BioModule')
import d00_sample as d00
import pysam

conn=mb.connect(host="localhost",user="root",passwd="123456",db='lihang')
cursor = conn.cursor()

conn2=mb.connect(host="162.105.59.51",port=8999,user="root",passwd="123456",db='hg19')
cursor2 = conn2.cursor()

#GET_RANDOM_1K_EXONS
#cursor.execute("SELECT chr,left_pos,right_pos FROM hg19_ensembl.exons WHERE id%100=1 and length(chr)<6")
#exons = cursor.fetchall()
cursor2.execute("SELECT chr,left,right FROM hg19.SureSelect_exome WHERE id%100=1 and length(chr)<6")
exons = cursor.fetchall()


#GET_COUNTS_OF_READS
samples = ['PC1_M','PC2_M','PC1_e','PC2_e','PT1_M','PT2_M','PT1_e','PT2_e']
samfiles = []
for sample in samples:
  file_path = d00.get_path('lihang',sample,'bam_bwa_genome_rmd')
  samfiles.append(pysam.Samfile(file_path, "rb"))
for region in exons:
  print region[0]+":"+str(region[1])+"-"+str(region[2]),region[2]-region[1],
  for samfile in samfiles:
    print samfile.count(region[0],region[1],region[2]),
  print ""
  
