import sys,time
sys.path.append('/home/fanxiaoying/lib/lib64/python/Bio')
sys.path.append('/data/Analysis/fanxiaoying/project/project01_polyA-RNAseq/modules')
import infra01_pos2info as in1
import MySQLdb as mb

def check_table_colume(cursor,conn,tablename,colume,info):
	cursor.execute("show columns from "+tablename)
	r0 = cursor.fetchall()
	lists = [i[0] for i in r0]
	if colume in lists:
		print "EXISTS"
	else:
		sql = "alter table %s add %s %s;" %(tablename,colume,info)
		cursor.execute(sql)
		conn.commit()

def append_colume_info_to_tables(cursor,conn,tablename,colume,info,ids,datas):
	check_table_colume(cursor,conn,tablename,colume,info)
	length_table = len(ids)
	values = [[datas[i],ids[i]] for i in range(length_table)]
	cursor.executemany("update "+tablename+" set "+colume+" = %s where id = %s ",values)
	conn.commit()
