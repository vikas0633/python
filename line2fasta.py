### line to fasta

import sys

count = 0
for line in open(sys.argv[1],'r'):
	line = line.strip()
	count += 1
	print '>'+str(count)+'_'+str(line)+'_'+str(len(line))
	print line