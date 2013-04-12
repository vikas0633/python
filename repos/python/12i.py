#12i.py - /plant/2011_week38/20110924/ for counting regulated sequences in cluster

import sys

	
### open cluster file
header=True
cluster_seq={}
for line in open(sys.argv[2],'r'):
	if(header==False):
		line=line.strip()
		token=line.split('\t')
		tokens=token[3].split(',')
		for i in range(1,len(tokens)):
			cluster_seq[tokens[i]]=1
	header=False

### open regulated sequence file
seq={}
count=0
header=True
for line in open(sys.argv[1],'r'):
	if(header==False):
		line=line.strip()
		token=line.split('\t')
		key=token[0]
		if key in cluster_seq:
			count += 1
		else:
			continue
	header=False
print "number of regulated sequences in the cluster:", (count)