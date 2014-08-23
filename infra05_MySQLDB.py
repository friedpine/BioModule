def CREATE_TABLE(cursor,tablename,colnames,types):
	print cursor,tablename,colnames,types
	sql = "CREATE TABLE "+tablename+" ("
	print sql
	for id,col in enumerate(colnames):
		print id,col
		sql += col+" "+types[id]+","
	sql = sql[0:len(sql)-1]+")"
	print sql
	try:
		cursor.execute(sql)
	except:	
		print "EXISTS"
	

