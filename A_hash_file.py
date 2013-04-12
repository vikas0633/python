###a script that hash the second column using first column as key

def hash_file(file):
	'''a script that hash the second column using first column as key'''
	hash = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 1:
			if line[0] != '#':
				token = line.split('\t')
				if len(token) > 1:
					hash[token[0]] = token[1]
	
	return hash
					