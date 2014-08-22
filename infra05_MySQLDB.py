def CREATE_TABLE(cursor,tablename,colnames,types):
  sql = "CREATE TABLE "+tablename+" ("
  for id,col in enumerate(colnames):
    sql += colnames+" "+types+","
  sql[-1] = ")"
  try:
    cursor.execute(sql)
  except:
    "EXISTS"

def GET_DATA_FROM_DB(cursor,conn,sql):
  return 1
