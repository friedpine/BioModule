class db()
	conn=mb.connect(host="localhost",user="root",passwd="123456",db=dbname)
	cursor = conn.cursor()
	def
	cursor.execute("select * from "+table1)