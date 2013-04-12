### 12l.py -/plant/2011_week42/20111017 - for making profile from fastq files

import sys

seq_count={} ### for counting read abundance
line_number=0 ### for tracking seq on 2nd line and so on
for line in open(sys.argv[1],'r'):
	line=line.strip()
	line_number += 1
	if( (line_number - 2)%4 == 0):
		if line in seq_count:
			seq_count[line] += 1
		else:
			seq_count[line] = 1
		
			
			
### open last profile file
seq = {}
for line in open(sys.argv[2],'r'):
	if(len(line) > 0):
		if(line[0] != '#'):
			line=line.strip()
			token=line.split('\t')
			seq[token[0]]=line ### read as key
			l=len(token)
		else:
			line=line.strip()
			token=line.split('\t')
			l=len(token)
			header=line

o=open(sys.argv[2],'w')
### be careful, it is only file specific
o.write(header+'\t'+str('lib-'+str(l))+'\n')
for key in seq_count:
	if key in seq:
		o.write(seq[key]+'\t'+str(seq_count[key])+'\n')
	else:	
		o.write(key+'\t')
		for i in range(1,l):
			o.write('0'+'\t')	
		o.write(str(seq_count[key])+'\n')
		
for key in seq:
	if key in seq_count:
		continue
	else:
		o.write(seq[key]+'\t'+'0'+'\n')