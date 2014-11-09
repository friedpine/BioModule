import MySQLdb as mb
import re,os

def table_2_dict(cursor,tablename,columes):
	gt = {}
	cursor.execute("select %s,%s from %s" %(columes[0],columes[1],tablename))
	r0 = cursor.fetchall()
	for i in r0:
		gt[i[0]] = i[1]
	return gt

def get_path2(cursor,sample,type):
	cursor.execute("select path from files where sample = %s and type = %s",([sample,type]))
	try:
		return cursor.fetchall()[0][0]
	except:
		return 'NA'

def get_path1(cursor,tablename,sample,filetype):
	sql = "select path from "+tablename+" where sample = '"+sample+"' and type = '"+filetype+"'"
	cursor.execute(sql)
	try:
		return cursor.fetchall()[0][0]
	except:
		return 'NA'

def get_colume_data(cursor,tablename,subjects,querys):
	sql = "select "+' '.join(subjects)+' from '+tablename+" where "+' and '.join(querys)
	print sql
	cursor.execute(sql)
	try:
		return cursor.fetchall()[0]
	except:
		print "NO SUCH THING!!!"
		return 'NA'

def get_ref(cursor,species,format,info,server="TANG"):
	try:
		cursor.execute("select path from bioinfo.ref where server=%s and species=%s and format=%s and info=%s",([server,species,format,info]))
		return cursor.fetchall()[0][0]
	except:
		print "Failed to get the REF files!"
		return 0

def get_script(cursor,type,server="TANG"): 
    try:
        cursor.execute("select path,cmds from bioinfo.scripts where server=%s and type=%s",([server,type]))
        return cursor.fetchall()[0]
    except: 
        print "Failed to get the SCRIPT files!"
        return 0

def check_validness_of_bamfiles(cursor,tablename,samples,filetype):
	bad_bam = []
	bad_index = []
	for sample in samples:
		bamfile = get_path1(cursor,tablename,sample,filetype)
		bamindex = bamfile+".bai"
		if not os.path.exists(bamfile):
			bad_bam.append(sample)
		if not os.path.exists(bamindex):
			bad_index.append(sample)
	if bad_index==[] and bad_bam==[]:
		print "Every_bam_and_index_files_exists!"
		return 1
	else:
		print "BAM_files_not_exists ",bad_bam
		print "INDEX_files_not_exists ",bad_index
		return 0
	
def samplename_transformer(cursor,conn,in_type,out_type,sample_name):
	try:
		cursor.execute("select "+out_type+" from samples where "+in_type+" = %s",([sample_name]))
		name = cursor.fetchall()[0][0]
	except:
		print "NO",out_type,"for",in_type,sample_name
		name = "NA_"+sample_name
	return name
	
def cp_move_files(pathes,operation,rec):
	if operation == 'cp':
		for i in pathes:
			subprocess.call("cp "+i[0]+" "+i[1],shell='True')
	if operation == 'mv':
		for i in pathes:
			subprocess.call("mv "+i[0]+" "+i[1],shell='True')	
	
def get_sample_file(cursor,sample,type):
	cursor.execute("select path from files where sample = %s and type = %s",[sample,type])
	return cursor.fetchall()[0][0]
	
def get_sample_info(cursor,sample,type):
	cursor.execute("select "+type+" from samples where sample = %s",[sample])
	return cursor.fetchall()[0][0]
	
def insert_sample_file(cursor,conn,sample,type,path):
	cursor.execute("insert ignore into files (sample,type,path,state)values(%s,%s,%s,NULL)",[sample,type,path])
	conn.commit()
