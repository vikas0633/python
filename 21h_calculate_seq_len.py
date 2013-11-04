# 21h_calculate_seq_len.py
# script simply takes fasta file as an input as prints
# out all the sequences name with their length in decreasing order



import os,sys, getopt
### function to calculate length of each seq in fasta file


def options(argv):
    infile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
       print 'python 21h_calculate_seq_len.py -i <inputfile>'
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print 'python 21h_calculate_seq_len.py -i <inputfile>'
          sys.exit()
       elif opt in ("-i", "--ifile"):
          infile = arg
    return infile
	
def fasta_len(file):
	first_line = True
	length = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0 :			
			if line[0] == '>':
				if first_line == False:
					length[header] = str_len
				str_len = 0
				header = line[1:]
			else:
				str_len += len(line)
		first_line = False			
	length[header] = str_len
	
	return length


def print_len(length):
	import operator
	sorted_len = sorted(length.iteritems(), key=operator.itemgetter(1),reverse=True)
	for key in sorted_len:
		print str(key[0]) +'\t'+ str(key[1])	


if __name__ == "__main__":
    
    infile = options(sys.argv[1:])    
    
    ### get the lengths back
    length = fasta_len(infile)
    
    ### print the values
    print_len(length)