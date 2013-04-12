## 21p_filter_fasta.py - script to filter fasta file based on the length of the sequences - /Users/vikas0633/Desktop/script/python

import os,sys, getopt
### function to calculate length of each seq in fasta file


def options(argv):
	infile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:l:",["ifile=","length="])
	except getopt.GetoptError:
		print 'python 21p_filter_fasta -i <inputfile> -l length'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'python 21p_filter_fasta -i <inputfile> -l length'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			infile = arg
		elif opt in ("-l", "--length"):
			length = int(arg.replace(',',''))
	return infile, length
	
def fasta_len(file):
	first_line = True
	seqs = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0 :
			if line[0] != '#':		
				if line[0] == '>':
					if first_line == False:
						seqs[header] = seq
					seq = ''
					header = line[1:]
				else:
					seq += line
				first_line = False			
	seqs[header] = seq
	
	return seqs


def print_seq(seqs,length):
	for key in sorted(seqs):
		if len(seqs[key]) >= length:
			print str('>'+key) 
			print str(seqs[key])	


if __name__ == "__main__":
    
    infile,length = options(sys.argv[1:])    
    
    ### get the lengths back
    seqs = fasta_len(infile)
    
    ### print the fasta file longer than length variable
    print_seq(seqs,length)