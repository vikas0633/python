### C_loadFasta.py - script to load fasta sequences 

def LOADfasta(file):
	first_line = True
	seq = {}
	string = ''
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0 :			
			if line[0] == '>':
				if first_line == False:
					if string != '': 
						seq[header] = string
				string = ''
				header = line[1:].strip().split()[0]
			else:
				string += line
		first_line = False			
	if string != '': 
		seq[header] = string
	return seq
