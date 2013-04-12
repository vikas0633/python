### scprit for couting correct hit in the positive control out put

import sys

count=0
hash={}
for line in open(sys.argv[1],'r'):
	line=line.strip()
	token=line.split('\t')
	if(int(token[8])<int(token[9])):
		st=int(token[8])
	else:
		st=int(token[9])
	###
	### hash the output file containing all top hits
	### take three elements as key (gap_pos, rep_name and insert_score)
	count += 1
	start=(10000000000*st)+count
	hash[int(start)]=token[0]
	
found=0	
for line in open(sys.argv[2],'r'):
	line=line.strip()
	token=line.split('\t')
	key=int(token[0])
	if key in hash:
		if(hash[key]==token[1]):
			found += 1
		else:
			print(line)
	
print (str(round((float(found)/float(count))*100,2))+'%')