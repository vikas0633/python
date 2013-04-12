### 12w.py filter out the sequences with minimum score cut-off

import sys


def filter_sequence(infile,cutoff):
	for line in open(infile,'r'):
		line=line.strip()
		token=line.split('\t')
		if(len(line)>0):
			if(line[0]=='#'): #header line
				i = token.index('score')
				print line
			else:
				if(float(token[i])>=float(sys.argv[2])):
					print line
if __name__ == "__main__":
	filter_sequence(sys.argv[1],int(sys.argv[2]))