### Usages the ORFfinder output and takes out sequences from a fasta file containing ORFs.
### perl orffinder.pl --infile=20121011_culinks.fa --outfile=20121011_culinks.orf.fa

### 21j_orf2fasta.py - script takes fasta file and output from orffinder and take out sequences with the orfs - /Users/vikas0633/Desktop/script/python


import os,sys, getopt

def options(argv):
    infile = ''; fasta_file = ''
    try:
    	opts, args = getopt.getopt(argv,"hi:f:",["orf_out=","fasta_file="])
    except getopt.GetoptError:
    	print 'python 21j_orf2fasta.py -i <orf_out> -f <fasta_file>'
    	sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print 'python 21j_orf2fasta.py -i <orf_out> -f <fasta_file>'
    		sys.exit()
    	elif opt in ("-i", "--orf_out"):
    		infile = arg
        elif opt in ("-f", "--fasta_file"):
        	fasta_file = arg
    return infile, fasta_file
    

def hash_orf(infile):
	orf = {}
	for line in open(infile,'r'):
		line = line.strip()
		token = line.split('\t')
		orf[token[0]]=''
	print len(orf)
	return orf
	
def print_fasta(orf,fasta_file):
	flag = False
	for line in open(fasta_file,'r'):
		line=line.strip()
		if line[0] == '>':
			if line[1:] in orf:
				flag = True
				orf.pop(line[1:])
			else:
				flag = False			
		if flag == True:
			print line
if __name__ == "__main__":
    
    infile,fasta_file = options(sys.argv[1:])    
    
    ### hash the orf-outfile
    orf = hash_orf(infile)
    
    ### print fasta
    print_fasta(orf,fasta_file)