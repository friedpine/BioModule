import re

def remove_low_quantity(seq,qua):
  pos = re.search("#*?$",qua).start()
  return seq[0:pos],qua[0:pos]
  
def reverse_complementary(seq):
  comp = {"A":"T","T":"A","C":"G","G":"C","N":"N"}

def adaptor_pair1_pos(seq):
  pos = [-1]
  try:
    pos.append(re.search("AGATCGGAAGA",seq).start())
    pos.append(re.search("GCACACGTCTG",seq).start()-11)
    pos.append(re.search("GGAAGAGCACA",seq).start()-5)
  except:
    tmp = 1
  return pos

def adaptor_pair2_pos(seq):
  pos = [-1]
  try:
    pos.append(re.search("AGATCGGAAGA",seq).start())
    pos.append(re.search("GCGTCGTGTAG",seq).start()-11)
    pos.append(re.search("GGAAGAGCGTC",seq).start()-5)
  except:
    tmp = 1
  return pos



def pairs_merge_by_adaptor(seq1,seq2,qua1,qua2):
  seq1,qua1 = remove_low_quantity(seq1,qua1)
  seq2,qua2 = remove_low_quantity(seq2,qua2)
  adp_pos_1 = adaptor_pair1_pos(seq1)
  adp_pos_2 = adaptor_pair2_pos(seq2)
  adp_pos = max(adp_pos_1+adp_pos_2)
  if adp_pos > 0:
    return 1,seq1[0:adp_pos],qua1[0:adp_pos]
  else:
    return 0,"NA","NA"

def pair_merge_by_overlap(seq1,seq2,qua1,qua2):
  

seq1 = "CACACCGGTCGCAGCTGCTCCGGCACCGGTACTGCCGGCGCGCCGGCAAACGTCTTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAAATCTC"
qua1 = "@@CFFFFFCCDDHGGHIIBHBGGIEHIGGIHHIIIIIHFCCBBBBBBCCCB@AAABACCCCB7?BB<AC>ACC>??9@CCCCCCCCC>CCBAAB0>:4:>C"
seq2 = "AAGACGTTTGCCGGCGCGCCGGCAGTACCGGTGCCGGAGCAGCTGCGACCGGTGTGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGG"
qua2 = "@@@?DDADFHFFF@@@DAG@DHIID@A>EEEDBACBBBB;9C?C@>9@<<BB5<6?3@CCC@B2?B@BCB@?B;B?CCCB?BB@A?ACCECAC>DA@B?BB"

seq1 = "CAAATCATATTATAGGCATGTGTGTGCAGTTGGTGGTGTATTAAATTTACTCACTATGTGGGCTGCAAACTTAGTAGGATTTGCTGTTGGGTTGGATGGTA"
seq2 = "AAAATCAAACTTATTTTATACTGACCATCTGACGTTCCAAAAATATTACTTAATAATGATTTCATACCATCCAACCCAACAGCAAATCCTACTAAGTTTGC"
pairs_merge_by_adaptor(seq1,seq2,qua1,qua2)
  



  
  
