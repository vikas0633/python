#21b_better_header.py - this script keeps only 4th field separated by '|'  - /Users/vikas0633/Desktop/plant/2012_week29/21b_better_header.py

import sys

file = sys.argv[1]

for line in open(file,'r'):
	line = line.strip()
	if line[0] == '>':
		token = line.split('|')
		anno = token[4].replace(' ','_')
		print token[0]+'_'+token[1]+'_'+anno
	else:
		print line
