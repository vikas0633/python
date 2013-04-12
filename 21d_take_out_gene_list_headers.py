
def fasta_len(file):
	first_line = True
	length = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0 :			
			if line[0] == '>':
				if first_line == False:
					length[header] = str_len
				str_len = ''
				header = line[1:]
			else:
				str_len += line
		first_line = False			
	length[header] = str_len
	
	return length

def extract_seq(file,header,length):
	for line in open(header,'r'):
		header = line.strip()	
		header = header.replace(' ','').replace('>','')
		if header in length:
			print '>'+header
			print length[header]
		
	
if __name__ == "__main__":
    
    # get the file and header
    file="/Users/vgupta/Desktop/temp/test.fa"
    header="/Users/vgupta/Desktop/temp/test.header"
    
    length = fasta_len(file)
    # get the sequence
    extract_seq(file,header,length)