### python 12k.py - /plant/2011_week40/20111004/ for counting unique size sequences and plotting the distribution from cluster files

import sys

seq={}
seq_len={}
for line in open(sys.argv[1],'r'):
	line=line.strip()
	if(len(line)>0):
		if(line[0]!='#'):
			token=line.split('\t')
			tokens=token[3].split(',')
			for i in range(1,len(tokens)):
				length=len(tokens[i])
				if tokens[i] in seq:
					seq[tokens[i]] += 1
				else:
					seq[tokens[i]] = 1
			for i in range(1,len(tokens)):
				length=len(tokens[i])
				if length in seq_len:
					seq_len[length] += 1
				else:
					seq_len[length] = 1
print "size fractionated reads count: " , (seq_len)	
print "total number of reads in cluster: " , (sum(seq.values()))
print "total number of unique reads in cluster: " , (len(seq))