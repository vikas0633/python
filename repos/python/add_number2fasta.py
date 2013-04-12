import sys

count = 0
for line in open(sys.argv[1],'r'):
	line = line.strip()
	
	if line[0]=='>':
		count += 1
		print line+'_'+str(count)
	else:
		print line


