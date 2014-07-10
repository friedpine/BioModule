import random as rd
def newfunc_set_count_sort(all):
	t = {}
	out = ''
	for i in all:
		if i not in t:
			t[i] = 1
		else:
			t[i] += 1
	n = sorted(t.items(),key=lambda d:d[1],reverse=True)
	for i in n:
		out = out+str(i[1])+'*'+i[0]+'#'
	return out
def ranges_regulate(range):
	if range[1]<range[0]:
		range = [range[1],range[0]]
	return range
def ranges_overlap(r1,r2):
	out = [0,0,0]
	r1 = [int(r1[0]),int(r1[1])]
	r2 = [int(r2[0]),int(r2[1])]
	if r2[0]<=r1[0] and r2[1]>=r1[0]:
		out[0] = 'T'
		map_len = min(r1[1],r2[1])-r1[0]
		out[1] = round(map_len/float(r1[1]-r1[0]),2)
		out[2] = round(map_len/float(r2[1]-r2[0]),2)
	elif r2[0]>=r1[0] and r2[0]<=r1[1]:
		out[0] = 'T'
		map_len = min(r1[1],r2[1])-r2[0]
		out[1] = round(map_len/float(r1[1]-r1[0]),2)
                out[2] = round(map_len/float(r2[1]-r2[0]),2)
	else:
		out = ['F',0,0]
	return out
def get_longest_range(ranges):
	a = sorted(ranges,key=lambda x:(x[1]-x[0]),reverse= True)	
	out = [a[0],ranges.index(a[0])]
	return out
def ranges_minus(big_range,little_ranges,least_gap):
	sorted_ranges = sorted(little_ranges)
	gap_ranges = [[sorted_ranges[i][1]+1,sorted_ranges[i+1][0]-1] for i in range(len(sorted_ranges)-1) if sorted_ranges[i][1]+least_gap<sorted_ranges[i+1][0]-1]
	if sorted_ranges[0][0]>big_range[0]:
		gap_ranges.append([big_range[0],sorted_ranges[0][0]-1])
	if sorted_ranges[-1][1]<big_range[1]:
		gap_ranges.append([sorted_ranges[-1][1]+1,big_range[1]])
	return sorted(gap_ranges)
def range_sort_len(ranges,reverse_or_not):
	a = sorted(ranges,key=lambda x:(x[1]-x[0]),reverse= reverse_or_not)
	return a
def random_sub_range(range,length):
	a = int(rd.random()*(range[1]-range[0]))
	return [a,a+length-1]
def merge_ranges(ranges,gap):
	ranges = sorted(ranges,key=lambda x:x[0])  
	for i in range(1,len(ranges)):
		if ranges[i][0]<=ranges[i-1][1]+gap:
			ranges[i] = [ranges[i-1][0],ranges[i][1]]
			ranges[i-1] = [0,0]
	return [i for i in ranges if i != [0,0]]