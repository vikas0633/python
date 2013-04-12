### for taking out pair where both reads mapped to same region
import sys

sam={} ###storing uniq lines

for line in open(sys.argv[1],'r'):
	line=line.strip()
	token=line.split('\t')
	key1=token[0].split('/')
	key2=token[2].split('_')
	key=key1[0],key2[0]
	if key in sam:
		sam[key] += 1
	else:
		sam[key] = 1

for line in open(sys.argv[1],'r'):
	line=line.strip()
	token=line.split('\t')
	key1=token[0].split('/')
	key2=token[2].split('_')
	key=key1[0],key2[0]

	if(sam[key] == 1):
		print line