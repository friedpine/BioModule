import sys,re,os,subprocess
import cPickle as pickle
import MySQLdb as mb
import time
db_name = sys.argv[1]
table_name = sys.argv[2]
bamfile = sys.argv[3]
dist = int(sys.argv[4])
record = sys.argv[5]
length_anchor = int(sys.argv[6])

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
			if pos_f-pos_s-length_anchor+1>0:
				out.append((record,chr,'pos',pos_s,pos_f+length_anchor-1,complete_seq,mate_seq,head))
		elif (flag & 0x10 != 0) and (flag & 0x20 != 0):
			if pos_s-pos_f-length_anchor+1>0:
				out.append((record,chr,'neg',pos_f,pos_s+length_anchor-1,complete_seq,mate_seq,head))
print record,len(out)

conn=mb.connect(host="localhost",user="root",passwd="123456",db=db_name)
cursor = conn.cursor()
sql = "CREATE TABLE "+ table_name+ """(
  `id` int(11) NOT NULL AUTO_INCREMENT,
  `sample` varchar(20) DEFAULT NULL,
  `chr` varchar(20) DEFAULT NULL,
  `strand` varchar(5) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `end` int(11) DEFAULT NULL,
  `seq` varchar(120) DEFAULT NULL,
  `mate` varchar(120) DEFAULT NULL,
  `head` varchar(100) DEFAULT NULL,
  `self_gencode` int(11) DEFAULT '0',
  `mate_gencode` int(11) DEFAULT '0',
  PRIMARY KEY (`id`),
  UNIQUE KEY `chr` (`chr`,`strand`,`start`,`end`,`sample`)
)
"""
try:
	cursor.execute(sql)
	conn.commit()
except:
	print "ALREADY"


cursor.executemany("insert ignore into "+table_name+" values(NULL,%s,%s,%s,%s,%s,%s,%s,%s,0,0) ",out)
conn.commit()
cursor.close()
conn.close()
