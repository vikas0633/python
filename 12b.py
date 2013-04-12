###for taking out sequences with particular pattern 

import sys, re

pattern=sys.argv[2]

for line in open(sys.argv[1],'r'):
	line=line.strip()
	if(line[0]==">"):
		seq_name=line
	if(line[0:3]==pattern):
		print seq_name
		print line
		