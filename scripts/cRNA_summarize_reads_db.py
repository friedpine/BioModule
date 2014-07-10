import sys,re,os,subprocess
import cPickle as pickle
import MySQLdb as mb
import time
db_name = sys.argv[1]
table_name = sys.argv[2]
bamfile = sys.argv[3]
dist = int(sys.argv[4])
record = sys.argv[5]

out = [] 
cmd = "samtools view -f 64 -F 12 "+bamfile
p1 = subprocess.Popen(cmd,shell = True,stdout=subprocess.PIPE)
for line in p1.stdout:
	array = str.split(line)
	flag = int(array[1])
	pos_f = int(array[3])
	pos_s = int(array[7])
	chr = array[2]
	if array[6] != '=' or abs(pos_s-pos_f)>dist:
		continue
	head = re.split(r'::',array[0])[0]
	complete_seq = re.split(r'::',array[0])[1]
	mate_seq = re.split(r'::',array[0])[2]
	if (array[6] == '=') and abs(pos_s-pos_f)<dist:
		if (flag & 0x10 == 0) and (flag & 0x20 == 0):
			if pos_f-pos_s-24>0:
				out.append((record,chr,'pos',pos_s,pos_f+24,complete_seq,mate_seq,head))
		elif (flag & 0x10 != 0) and (flag & 0x20 != 0):
			if pos_s-pos_f-24>0:
				out.append((record,chr,'neg',pos_f,pos_s+24,complete_seq,mate_seq,head))
print record,len(out)

conn=mb.connect(host="localhost",user="root",passwd="123456",db=db_name)
cursor = conn.cursor()
# for i in range(int(len(out)/10)+1):
# 	print record,i,'B1',time.time()
# 	cursor.executemany("insert ignore into "+table_name+" values(NULL,%s,%s,%s,%s,%s,%s,%s,%s) ",out[100000*i:min(100000*i+100000,len(out))])
# 	print record,i,'B2',time.time()
# 	conn.commit()
# 	print record,i,'B3',time.time()
i = 1
print record,i,'B1',time.time()
cursor.executemany("insert ignore into "+table_name+" values(NULL,%s,%s,%s,%s,%s,%s,%s,%s,0,0) ",out)
print record,i,'B2',time.time()
conn.commit()
print record,i,'B3',time.time()
cursor.close()
conn.close()
