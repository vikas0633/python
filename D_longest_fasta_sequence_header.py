#D_longest_fasta_sequence_header.py - /Users/vikas0633/Desktop/script/python/ - script return headers of longest sequence

def longest_seq(file):
	first_line = True
	hash = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0: 
			if line[0] != '#':
				if line[0] == '>': ## new sequence
					
					### hash the length
					if (first_line == False):
						hash[header] = seq_len
					header = line
					seq_len = 0
					token = line.split('.')
					tokens = line.split(':')
					seq_strand = tokens[4].strip()
					first_line = False
				else:
					seq_len += len(line)

	
	### for last sequence
	hash[header] = seq_len
	
	
	### find the longest length
	if len(hash) == 0:
		hash['fake'] = ' '
	return max(hash, key=hash.get)