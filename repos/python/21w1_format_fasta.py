

### check output with 
# python ~/script/python/21w1_format_fasta.py -i Ljr_cdna.chr0-6.fa.refined -g 

import os,sys,getopt, re
from C_loadFasta import *

### main argument to 

def options(argv):
	inputfile = ''
	gff3 = ''
	try:
		opts, args = getopt.getopt(argv,"hi:g:",["ifile=","gff3="])
	except getopt.GetoptError:
		print 'python 21w1_format_fasta.py -i <inputfile> -g <gff3>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'python 21w1_format_fasta.py -i <inputfile> -g <gff3>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-g", "--gff3"):
			gff3 = arg
	return inputfile, gff3
    
def HASHfasta(inf):
	return LOADfasta(inf)
	
def printFasta(gff3,hash):
	uniq = {} ### to write only unique sequences
	for line in open(gff3,'r'):
		line = line.strip()
		token = line.split('\t')
		if len(line) > 1:
			if line[0] != '#':
				if token[2] == "mRNA":
					if token[1]=="CUFFLINKS":
						match = re.search(r'ID=.+;',line)
						match = match.group().split(';')[0].replace('ID=','')
						if match not in uniq:
							if match in hash:
								uniq[match] = ''
								print '>'+match
								print hash[match]
					else:
						match = re.search(r'ID=.*',line)
						match = match.group().split(';')[0].replace('ID=','')
						#match = match.group().split(';')[0].replace('Parent=','')
						if match not in uniq:
							if match in hash:
								uniq[match] = ''
								print '>'+match
								print hash[match]
							
		

if __name__ == "__main__":
    
	inf,gff3 = options(sys.argv[1:])
	
	hash = HASHfasta(inf)
	
	### print duplicate removed fasta
	printFasta(gff3,hash)