###21d_take_out_gene.py - this script takes out a sequence from fasta file given correct header name 
### nice -n 19 python 21d_take_out_gene.py -f <fasta file> -n <header_name>
### python 21d_take_out_gene.py -f Kasuza.fasta.refined -n chr5_31999127_31999213

import sys, getopt

def file_empty(file):
    count = sum([1 for line in open(file)])
    if count == 0:
        sys.exit(file+' is empty')

### get the options 
def options(argv):
    try:
       opts, args = getopt.getopt(argv,"hf:n:",["fasta=","header="])
    except getopt.GetoptError:
       print 'python 21d_take_out_gene.py -f <fasta file> -n <header_name> '
       sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print 'python 21d_take_out_gene.py -f <fasta file> -n <header_name>'
    		sys.exit()
    	elif opt in ("-f", "--fasta"):
    		file = arg
    	elif opt in ("-n", "--header"):
       		header = arg.replace('"','')
    return file, header

def extract_seq(file,header):
	header = header.replace(' ','').replace('>','')
	flag = False
	for line in open(file,'r'):
		if line[0]=='>':
			line = line.strip().split(' ')[0]
			if line[1:].strip() == header:
				print line
				flag = True
			else:
				flag = False
		else:
			if flag == True:
				line = line.strip()
				print line
		
	
if __name__ == "__main__":
    
    # get the file and header
    file,header = options(sys.argv[1:])
    
    file_empty(file)
    # get the sequence
    extract_seq(file,header)
        
    
