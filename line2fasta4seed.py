### line to fasta

import sys

count = 0

def add_N(x):
	string = ''
	for i in range(x):
		string += 'A'
	return string
for line in open(sys.argv[1],'r'):
	line = line.strip()
	count += 1
	if len(line) < 21:
		x = 21 - len(line) 
		line += add_N(x)
	print '>'+str(count)
	print line