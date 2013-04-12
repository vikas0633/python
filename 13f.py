### for controling signle end directionality

import sys


nb=int(sys.argv[3]) ### number of one side hanging bases
mate_pair=sys.argv[4] ### insert library
### first read in pair
if(mate_pair=='F'):
	for line in open(sys.argv[1],'r'):
		line=line.strip()
		token=line.split('\t')
		tokens=token[2].split(',')
		up_stream=int(tokens[1])
		gap=int(tokens[3])
		if(int(token[1])==0): ### forward strand of query
			if(int(token[3]) <= up_stream):
				print(line)
		else:	### reverse strand of query
			if(int(token[3]) > up_stream):
				print(line)
		
		
	for line in open(sys.argv[2],'r'):
		line=line.strip()
		token=line.split('\t')
		if(int(token[1])==0): ###  forward strand of query
			if(int(token[3]) <= up_stream):
				print(line)
		else:	### reverse strand of query
			if(int(token[3]) > up_stream):
				print(line)
				
if(mate_pair=='T'):
	for line in open(sys.argv[1],'r'):
		line=line.strip()
		token=line.split('\t')
		print(line)
	for line in open(sys.argv[2],'r'):
		line=line.strip()
		token=line.split('\t')
		print(line)
