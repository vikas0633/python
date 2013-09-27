
def fasta_len(file):
	first_line = True
	seq = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0 and not line.startswith('#'):			
			if line[0] == '>':
				if first_line == False:
					seq[header] = sq
				sq = ''
				header = line[1:].split(' ')[0].split('|')[2]
			else:
				sq += line
			first_line = False			
	seq[header] = sq
	
	return seq

def extract_seq(file,header,seq):
	for line in open(header,'r'):
		header = line.strip()	
		header = header.replace(' ','').replace('>','')
		if header in seq:
			print '>'+header
			print seq[header]
		
	
if __name__ == "__main__":
    
    # get the file and header
    file="/Users/vgupta/Desktop/03_Lotus_annotation/Gma.seq.all"
    header="/Users/vgupta/Desktop/03_Lotus_annotation/svend/4_2.txt"
    
    seq = fasta_len(file)
    # get the sequence
    extract_seq(file,header,seq)