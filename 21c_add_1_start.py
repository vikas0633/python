#21c_add_1_start.py - this script can add +1 to start position in a fasta file
# python 21c_add_1_start.py <infile> > <outfile>
# should work with fasta file by bedtools

import sys

def add_1(file):
	for line in open(file,'r'):
		line = line.strip()
		if line[0] == '>':
			token=line.replace(':',' ').replace('-',' ').split()
			print token[0]+'_'+str(int(token[1])+1)+'_'+token[2]
		else:
			print line



if __name__ == "__main__":
	add_1(sys.argv[1])