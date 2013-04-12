#13q.py - 	take out all the contigs alraeady placed in the psuedomolecule

import sys 

### contig containing column
column=int(sys.argv[2])

### store contigs in hash
contig={}


### open file with placed contigs
for line in open(sys.argv[1],'r'):
	if(len(line)>0):
		line=line.strip()
		token=line.split('\t')
		if(len(token[column-1])>0):
			contig[token[column-1]]=1
			
### open contig file and print only once which are not in contig hash table
found=False
for line in open(sys.argv[3],'r'):
	line=line.strip()
	if(len(line)>0):
		if(line[0]=='>'):
			if line[1:] in contig:
				found=False
			else:
				found=True
				print line
				
		if((line[0]!='>') & (found==True)):
			print line
				
