import re

def remove_low_quantity(seq,qua):
  pos = re.search("#*?$",qua).start()
  return seq[0:pos]
def Adaptor_infos(seq,adaptor):
  
  
  
