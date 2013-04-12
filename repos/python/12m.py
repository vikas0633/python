### Distinguish between mapped and unmapped reads and return respective .f*q files
#Usage: python 12m.py mapped.sam infile.fq *_repeat_filtered > *_repeat 

import sys

### hash the sam file
seq_id={}
for line in open(sys.argv[1],'r'):
	line=line.strip()
	token=line.split('\t')
	seq_id[token[0]]=1
	

### open fastq file and output non-mapped reads
o=open(sys.argv[3],'w')


key_found=False
line_no=0
for line in open(sys.argv[2],'r'):
	line=line.strip()
	line_no += 1
	if((key_found==True) & (line[0]=='@') & (line_no%4==1)):
		key_found=False ### terminate printing
	if((line[0]=='@')& (line_no%4==1)): ### read id
		key=line[1:]
		if key in seq_id: ### check if id exist in the sam file hash
			key_found=True
			print line
		else:
			o.write(line);o.write('\n')
	if((key_found==True) & (line_no%4!=1)):
		print line
	if((key_found==False) & (line_no%4!=1)):
		o.write(line);o.write('\n')

o.close()