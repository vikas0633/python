##13n.py - 	correct sequence name in fasta file (remove every thing after spaces), problem when mapping

import sys

for line in open(sys.argv[1],'r'):
	line=line.strip()
	if(len(line) > 0):
		if(line[0]=='>'):
			token=line.split(' ')
			print (token[0])
		else:
			print (line)