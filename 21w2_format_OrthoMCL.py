#21w1_format_orthoMCL.py - format OrthoMCL output - /Users/vikas0633/Desktop/script/python
# idea is to make a column consisting unique mRNAids and corresponding OrthoMCL group in the next columns


# Usage: python /Users/vikas0633/Desktop/script/python/21w2_format_orthoMCL.py -i mclGroups.txt

### check output with 
# python /Users/vikas0633/Desktop/script/python/21w1_format_fasta.py -i Ljr_cdna.chr0-6.fa.refined| sort | uniq -c | sort -nr |head

import os,sys,getopt, re


### main argument to 

def options(argv):
    inputfile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:",["ifile=",])
    except getopt.GetoptError:
       print 'python 21w1_format_orthoMCL.py -i <inputfile> '
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print 'python 21w1_format_orthoMCL.py -i <inputfile> '
          sys.exit()
       elif opt in ("-i", "--ifile"):
          inputfile = arg
    
    return inputfile
    
def format_orthoMCL(inf):
	hash = {}
	for line in open(inf,'r'):
		line = line.strip()
		tokens = line.split(' ')
		for token in tokens:
			key = token.split('|')[0]
			if key == 'ljr':
				print token[4:]+'\t'+tokens[0][:-1]
			
	return hash
	
def printFasta(hash):
	for key in hash:
		print key
		print hash[key]


if __name__ == "__main__":
    
    inf = options(sys.argv[1:])
    
    hash = format_orthoMCL(inf)
    
    ### print duplicate removed fasta
    printFasta(hash)