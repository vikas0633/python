### script for combining sequences fro profile with the mapped tasiRNA sequences 

import sys

profile=sys.argv[1]
tasiRNA=sys.argv[2]
o=open(sys.argv[3],'w')

### load miRNA file
hash={}
pos={} ### to make sure that only one position is take and removes redundancy which comes from ta-siRNA standard output
first_line=True
for line in open(tasiRNA,'r'):
	line = line.strip()
	if(len(line)>0):
		# header
		if (first_line==True):
			first_line=False
		else:
			token = line.split('\t') 
			if token[3] in hash:
				if (token[0], token[1]) not in pos:
					pos[token[0],token[1]]=''
					hash[token[3]] += (','+token[0]+'_'+token[1])
			else:
				pos[token[0],token[1]]=''
				hash[token[3]] = 'predicted as ta-siRNA, position: '+token[0]+'_'+token[1]
### go through each sequence in profiles
first_line=True
for line in open(profile,'r'):
 	line =line.strip()
 	token=line.split('\t')
	if (first_line==True):
		first_line=False
	else:
		if token[0] in hash:
			o.write(token[0]+'\t'+hash[token[0]]+'\n')
		else:
			o.write(token[0]+'\t'+'not predicted as ta-siRNA'+'\n')
 		
o.close()