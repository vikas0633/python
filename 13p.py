#13p.py - 	script to remove pair mapped 

import sys

sam1={}

for line in open(sys.argv[1],'r'):
	line = line.strip()
	token=line.split('\t')
	id=token[0].split('/')
	seq_id=id[0]
	### store ids in SAM1
	sam1[seq_id]=1



sam_both={}
### open second sam file and look for best hit
for line in open(sys.argv[2],'r'):
	line = line.strip()
	token=line.split('\t')
	id=token[0].split('/')
	key=id[0]
	### store ids in SAM1
	if key in sam1:
		sam_both[key]=1
	else:
		continue
		
### open first sam and print output
o1=open(sys.argv[3],'w')
for line in open(sys.argv[1],'r'):
	line = line.strip()
	token=line.split('\t')
	id=token[0].split('/')
	key=id[0]
	if key in sam_both:
		continue
	else:
		o1.write(line)
		o1.write('\n')
		
### open second sam and print output
o2=open(sys.argv[4],'w')
for line in open(sys.argv[2],'r'):
	line = line.strip()
	token=line.split('\t')
	id=token[0].split('/')
	key=id[0]
	if key in sam_both:
		continue
	else:
		o2.write(line)
		o2.write('\n')