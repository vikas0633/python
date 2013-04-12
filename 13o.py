### script to take out a sequence from a multi fasta file
### Usage: nice -n 19 python 13o.py contig_file contig _name >contig.txt

import sys


seq_name=sys.argv[2]
line_found=False

for line in open(sys.argv[1],'r'):
	line=line.strip()
	if(len(line)>0):
		if(line[0]=='>'):
			if(seq_name==line[1:]):
				line_found=True
				print line
			else:
				line_found=False
		if(line[0]!='>'):
			if(line_found==True):
				print line


