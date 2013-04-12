### script for combining sequences fro profile with the mapped genomic sequences which contain genomic annotation

import sys

profile=sys.argv[1]
region=sys.argv[2]
o=open(sys.argv[3],'w')

### load region file in a hash using sequence as key
hash={}

first_line=True
for line in open(region,'r'):
	line = line.strip()
	# header
	if (first_line==True):
		first_line=False
	else:
		token = line.split('\t') 
		if token[0] in hash:
			hash[token[0]] = hash[token[0]] + ',' + token[1]
		else:
			hash[token[0]] = token[1]


### go through each sequence in profiles
for line in open(profile,'r'):
 	line =line.strip()
 	
 	token =line.split('\t')
 	
 	if token[0] in hash:
 		o.write(token[0]+'\t'+hash[token[0]]+'\n')
 	else:
 		o.write(token[0]+'\t'+'not mapped'+'\n')
 		
o.close()