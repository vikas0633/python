### script for combining sequences fro profile with the mapped miRNA sequences 

import sys

profile=sys.argv[1]
miRNA=sys.argv[2]
o=open(sys.argv[3],'w')

### load miRNA file
hash={}

first_line=True
for line in open(miRNA,'r'):
	line = line.strip()
	if(len(line)>0):
		# header
		if (first_line==True):
			first_line=False
		else:
			token = line.split('\t') 
			if token[6] in hash:
				hash[token[6]] = 'predicted as miRNA'
			else:
				hash[token[6]] = 'predicted as miRNA'


### go through each sequence in profiles
for line in open(profile,'r'):
 	line =line.strip()
 	
 	token =line.split('\t')
 	
 	if token[0] in hash:
 		o.write(token[0]+'\t'+hash[token[0]]+'\n')
 	else:
 		o.write(token[0]+'\t'+'not predicted as miRNA'+'\n')
 		
o.close()