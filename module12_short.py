import re

def remove_low_quantity(seq,qua):
  pos = re.search("#*?$",qua).start()
  return seq[0:pos],qua[0:pos]
  
def reverse_complementary(seq):
  comp = {"A":"T","T":"A","C":"G","G":"C","N":"N"}
  RC_seq = ""
  for i in seq:
    RC_seq += comp[i]
  return RC_seq[::-1]
  
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

def seq_overlaption(seq1,seq2,qua1,qua2):
  print "NICE!"

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
  seq1,qua1 = remove_low_quantity(seq1,qua1)
  seq2,qua2 = remove_low_quantity(seq2,qua2)  
  seq2 = reverse_complementary(seq2)
  qua2 = qua2[::-1]
  poses = range(0,len(seq2)-12,6)
  segs = [seq2[i:12+i] for i in poses]
  overlap_pos = [-1]
  for id,seg in enumerate(segs):
    try:
      overlap_pos.append(re.search(seg,seq1).start()-poses[id])
    except:
      tmp = 1
  pos = max(overlap_pos)
  if pos >= 0:
    return 1,seq1+seq2[len(seq1)-pos:len(seq2)],qua1+qua2[len(seq1)-pos:len(seq2)]
  else:
    return 0,seq1,seq2
  
seq1 = "CACACCGGTCGCAGCTGCTCCGGCACCGGTACTGCCGGCGCGCCGGCAAACGTCTTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAAATCTC"
qua1 = "@@CFFFFFCCDDHGGHIIBHBGGIEHIGGIHHIIIIIHFCCBBBBBBCCCB@AAABACCCCB7?BB<AC>ACC>??9@CCCCCCCCC>CCBAAB0>:4:>C"
seq2 = "AAGACGTTTGCCGGCGCGCCGGCAGTACCGGTGCCGGAGCAGCTGCGACCGGTGTGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGG"
qua2 = "@@@?DDADFHFFF@@@DAG@DHIID@A>EEEDBACBBBB;9C?C@>9@<<BB5<6?3@CCC@B2?B@BCB@?B;B?CCCB?BB@A?ACCECAC>DA@B?BB"

seq1 = "CAAATCATATTATAGGCATGTGTGTGCAGTTGGTGGTGTATTAAATTTACTCACTATGTGGGCTGCAAACTTAGTAGGATTTGCTGTTGGGTTGGATGGTA"
seq2 = "AAAATCAAACTTATTTTATACTGACCATCTGACGTTCCAAAAATATTACTTAATAATGATTTCATACCATCCAACCCAACAGCAAATCCTACTAAGTTTGC"

#pair_merge_by_overlap(seq1,seq2,qua1,qua2)
#pairs_merge_by_adaptor(seq1,seq2,qua1,qua2)
  



  
  
