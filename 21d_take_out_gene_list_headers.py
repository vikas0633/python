import re, sys

def fasta_len(file):
	first_line = True
	seq = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0 :			
			if line[0] == '>':
				if first_line == False:
					seq[header] = sq
				sq = ''
				header = line[1:].split()[0]
			else:
				sq += line
		first_line = False			
	seq[header] = sq
	return seq

def extract_seq(file,header_file,seq):
	o=open(header_file+'.fa','w')
	count = 0
	for line in open(header_file,'r'):
		line = line.strip()
		if len(line) > 0 and not line.startswith('#'):
			header = line.strip().split()[0]	
			header = header.replace(' ','').replace('>','')
			count += 1
			if header in seq:
				print str(count)+'. '+header
				o.write('>'+header+'\n')
				o.write(seq[header]+'\n')
				
			else:
				print 'Warning: Sequence not found for: ', header
		
	
if __name__ == "__main__":
    
    # get the file and header
    file=sys.argv[1]
    header=sys.argv[2]
    
    seq = fasta_len(file)
    # get the sequence
    extract_seq(file,header,seq)