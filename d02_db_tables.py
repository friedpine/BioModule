import sys,time
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import infra01_pos2info as in1
import MySQLdb as mb
import re


def check_table(cursor,database,tablename):
	pass

def check_table_columes(cursor,tablename,columes):
	try:
		cursor.execute("show columns from "+tablename)
		r0 = cursor.fetchall()
		lists = [i[0] for i in r0]
		out = []
		for colume in columes:
			if colume not in lists:
				out.append(colume)
		if out == []:
			print "The table is already Exists with every Colume!"
			return 1
		else:
			print "The table is already Exists but lacks some Columes!",out
			return 2
	except:
		print "TABLE DOES NOT EXISTS!"
		return 0

def add_table_colume(cursor,conn,tablename,colume,info):
	cursor.execute("show columns from "+tablename)
	r0 = cursor.fetchall()
	lists = [i[0] for i in r0]
	if colume in lists:
		print "Colume_Already_EXISTS!"
	else:
		sql = "alter table %s add %s %s;" %(tablename,colume,info)
		cursor.execute(sql)
		conn.commit()

def append_colume_info_to_tables(cursor,conn,tablename,colume,info,ids,datas):
	add_table_colume(cursor,conn,tablename,colume,info)
	length_table = len(ids)
	values = [[datas[i],ids[i]] for i in range(length_table)]
	cursor.executemany("update "+tablename+" set "+colume+" = %s where id = %s ",values)
	conn.commit()

def Create_Mysql_tables(cursor,conn,tablename,columes,infos):
	if len(columes) != len(infos):
		print "The input infos: columes and infos are not of the same length!"
		return
	sqlpara = ','.join([columes[x]+" "+infos[x] for x in range(len(columes))])
	sqlcmd = "create table %s(%s)" %(tablename,sqlpara)
	checkstate = check_table_columes(cursor,tablename,columes)
	if checkstate==1:
		return 1
	elif checkstate==0:
		try:
			cursor.execute(sqlcmd)
			return 1
		except:
			print "THE sqlcommand generate Errors while Executing!",sqlcmd 
			return 0 
	elif checkstate==2:
		print "The Table need altered!"
		return 0

def Insert_lines_Into_table(cursor,conn,tablename,columes,data):
	if data == []:
		print "NO_DATA_PROVIDED!"
		return 1 
	if len(data[0]) != len(columes):
		print "Columes and data not of the same length!"
		return 0
	columes_joined = ",".join(columes)
	percent_s = ",".join(["%s"]*len(columes))
	sqlcmd = "replace into %s (%s)values(%s)" %(tablename,columes_joined,percent_s)
	print sqlcmd
	cursor.executemany(sqlcmd,data)
	conn.commit()

		
def Save_text_to_Mysql(cursor,conn,textpath,tablename,columes,infos):
	table_state = Create_Mysql_tables(cursor,conn,tablename,columes,infos)
	if table_state == 0:
		return
	try:
		file = open(textpath,'r')
	except:
		print "Errors:File Reading"
		return
	values = []
	count = 0
	for line in file:
		count += 1
		array = re.split('\s+',line)
		values.append(array[:-1])
		if count%1000 == 0:
			Insert_lines_Into_table(cursor,conn,tablename,columes,values)
			print count
			values = []
	Insert_lines_Into_table(cursor,conn,tablename,columes,values)

def Append_text_to_Table(cursor,conn,textpath,tablename,columes):
	table_state = check_table_columes(cursor,tablename,columes)	
	if table_state == 1:
		try:
			file = open(textpath,'r')
		except:
			print "Errors:File Reading"
			return
		values = []
		for line in file:
			array = re.split('\s+',line)
			values.append(array[:-1])
		Insert_lines_Into_table(cursor,conn,tablename,columes,values)
		print "APPEND",len(values),"Lines into the table!"
