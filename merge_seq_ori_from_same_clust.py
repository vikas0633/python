### merge sequences originating from same clusters

import sys

hash_pos = {}

for line in open(sys.argv[1],'r'):
	if line[0] != '#':
		found = False
		line = line.strip()
		line = line.replace('"','')
		tokens = line.split('\t')[1].split(';')
		if len(tokens) > 1:  
			token = tokens[1].split(',')
			for pos in token:
				key = pos.split('_')[0],pos.split('_')[1]
				if (key[0],int(key[1])) in hash_pos:
					found = True
					break
				else:
					for i in range(-20,20):
						hash_pos[key[0],int(key[1])+i] = ''
					found = False
		print line + "\tMatched Previous Loci: "+str(found)