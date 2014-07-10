import re, sys, os, copy
import subprocess
import time
import cPickle as pickle
import MySQLdb as mb
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import infra00_ranges_operate as in0


DB_NAME = sys.argv[1]
#table_exons = sys.argv[2]
# table_t_ex = sys.argv[3]
# table_transc = sys.argv[4]
# table_ref = sys.argv[5]


conn=mb.connect(host="localhost",user="root",passwd="123456",db=DB_NAME)
cursor = conn.cursor()
cursor.execute("SELECT id,mate FROM 11_raw WHERE id IN (SELECT id FROM 11_mm_exon where t_dis=0)")
for i in cursor.fetchall():
        print ">"+str(i[0])
        print i[1]


# exons = []
# t_ex = []
# transc = []

if "STEP1" in sys.argv:
        psl = sys.argv[1]
        file = open(psl)
        for i in range(5):
                file.readline()
        for line in file:
                t = re.split('\s+',line)
                if len(t)<10 or int(t[0])<0.95*int(t[10]):
                        continue
                chr = t[13]
                exon_count = int(t[17])
                exons_start = [int(i) for i in re.split(',',t[20])[:-1]]
                exons_size = [int(i) for i in re.split(',',t[18])[:-1]]
                exons_ranges = [[exons_start[i],exons_start[i]+exons_size[i]-1] for i in range(exon_count)]
                exons_ranges = in0.merge_ranges(exons_ranges,25)
                info = ''
                exonsize = ''
                intronsize = ''
                for i in exons_ranges:
                        exons.append([chr,i[0],i[1],i[1]-i[0]])
                        t_ex.append([t[9],chr,i[0],i[1]])
                        info += str(i[0])+'-'+str(i[1])+','
                        exonsize += str(i[1]-i[0])+','
                for i in range(1,len(exons_ranges)):
                        intronsize += str(exons_ranges[i][0]-exons_ranges[i-1][1]-1)+','
                transc.append([t[9],t[10],t[13],t[15],t[16],len(exons_ranges),info,exonsize,intronsize])
        print transc[-1]
        # cursor.executemany("insert ignore into "+table_exons+" (chr,left_pos,right_pos,length)values(%s,%s,%s,%s)",exons)
        # cursor.executemany("insert ignore into "+table_t_ex+" (transc,chr,left_pos,right_pos)values(%s,%s,%s,%s)",t_ex)
        # cursor.executemany("insert ignore into "+table_transc+" (id,length,chr,left_pos,right_pos,exon_count,info,exon_size,intron_size)values(%s,%s,%s,%s,%s,%s,%s,%s,%s)",transc)
        # conn.commit()