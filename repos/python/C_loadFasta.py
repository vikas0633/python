### C_loadFasta.py - script to load fasta sequences 

def LOADfasta(file):
	first_line = True
	seq = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0 :			
			if line[0] == '>':
				if first_line == False:
					if str != '': 
						seq[header] = str
				str = ''
				header = line[1:].strip().lstrip()
			else:
				str += line
		first_line = False			
	if str != '': 
		seq[header] = str
	return seq