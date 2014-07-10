import sys,re,os,subprocess
import cPickle as pickle
outdat = sys.argv[1]
bamfile = sys.argv[2]
dist = int(sys.argv[3])
record = sys.argv[4]

out = {}
out[record] = {}
out['all_pos_ff'] = 0
out['all_pos_fs'] = 0
out['all_neg_ff']= 0
out['all_neg_fs'] = 0
out['one_pos_neg'] = 0
out['duplication'] = 0 
out['total'] = 0

cmd = "samtools view -f 64 -F 12 "+bamfile
p1 = subprocess.Popen(cmd,shell = True,stdout=subprocess.PIPE)
for line in p1.stdout:
	array = str.split(line)
	flag = int(array[1])
	pos_f = int(array[3])
	pos_s = int(array[7])
	if array[6] != '=' or abs(pos_s-pos_f)>dist:
		continue
	head = re.split(r'::',array[0])[0]
	complete_seq = re.split(r'::',array[0])[1]
	matepos = re.split(r'::',array[0])[2]
	if (array[6] == '=') and abs(pos_s-pos_f)<dist:
		out['total'] += 1
		if array[2] not in out[record]:
			out[record][array[2]] = {}
		if (flag & 0x10 == 0) and (flag & 0x20 == 0):
			if pos_s-pos_f-24>0:
				out['all_pos_ff'] += 1
			elif pos_f-pos_s-24>0:
				out['all_pos_fs'] += 1
				if pos_s in out[record][array[2]]:
					out['duplication'] += 1
				out[record][array[2]][pos_s] = {}
				out[record][array[2]][pos_s]['end'] = pos_f+24
				out[record][array[2]][pos_s]['strand'] = 'pos'
				out[record][array[2]][pos_s]['seq'] = complete_seq
				out[record][array[2]][pos_s]['head'] = head
				out[record][array[2]][pos_s]['event'] = ''
				out[record][array[2]][pos_s]['mate'] = matepos
		elif (flag & 0x10 != 0) and (flag & 0x20 != 0):
			if pos_s-pos_f-24>0:
				out['all_neg_ff'] += 1
				if pos_f in out[record][array[2]]:
					out['duplication'] +=  1
				out[record][array[2]][pos_f] = {}
				out[record][array[2]][pos_f]['end'] = pos_s+24
				out[record][array[2]][pos_f]['strand'] = 'neg'
				out[record][array[2]][pos_f]['seq'] = complete_seq
				out[record][array[2]][pos_f]['head'] = head
				out[record][array[2]][pos_f]['event'] = ''
				out[record][array[2]][pos_f]['mate'] = matepos
			elif pos_f-pos_s-24>0:
				out['all_neg_fs'] += 1
		else:
			out['one_pos_neg'] += 1
pickle.dump(out,open(outdat,'w'))
