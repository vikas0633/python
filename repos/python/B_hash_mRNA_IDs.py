### returns a uniq mRNA id hash

import re
def hash_mRNA(file):
	ID = {}
	for line in open(file,'r'):
		token = line.strip().split()
		if token[2]== 'gene':
			### get the mRNA ID
			match = re.search(r'ID=.+;',line)
			if match:
				match = match.group().split(';')[0].replace('ID=','')
				ID[match]=''	
	return ID