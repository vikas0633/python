### combining results


import sys

    
#make mismatch file as hash to store values of read in hash
hash={}
for line in open(sys.argv[2],'r'):
	line=line.strip()
	token=line.split('_')
	try:
		hash[token[1]]=line
	except:
		continue
	
header=True
line_number=0
for line in open(sys.argv[1],'r'):
	line=line.strip()
	token=line.split()
	line_number += 1
	key=token[0] ### check if key in hash
	if(header==True):
		print str(line+'\t'+str('line_number')+'\t'+str('hits with mismatches(0:1:2:3)')+'\t'+'best_hit'+'\t'+'best_family')
	else:
		if key in hash:
			print str(line+'\t'+hash[key])
		else:
			print str(line+'\t'+str(line_number)+'\t'+str('0:0:0:0')+'\t'+'none'+'\t'+'none')
	header=False