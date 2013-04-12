### remove U and replace by T

import sys

for line in open(sys.argv[1],'r'):
	line = line.strip()
	if(line[0]=='>'):
		print(line)
	else:
		line=line.replace("U","T")
		print(line)