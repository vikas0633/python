# 10a.py script for reading genome file and taking out N regions one by one

import sys
import numpy as np
# open genome file

seq = '' ; a = np.array([])
pos = 0; last_pos = 0; len_rep = 0
gap_found=False
for line in open(sys.argv[1],'r'):
	line = line.strip()
	if(line[0] =='>'):
		seq = '' ; a = np.array([])
		pos = 0; last_pos = 0; len_rep = 0
		chr = line[1:]
	else:
		#split line if there is any N
		token = line.split('NNNN')
		if((len(token) > 1) &(pos != last_pos)):
			len_rep += (len(line) - len(token[0]))
			seq += line[len(token[0]):len(line)]
			pos += len(token[0])
			gap_found=True
			last_pos = pos
		elif(gap_found==True):
			if((len(token[len(token)-1]) >1) &(line[len(line)-1] != 'N')):
				len_rep += (len(line) - len(token[len(token)-1]))
				seq += line[0:(len(line) - len(token[len(token)-1]))]
				#print("nice -n 19 fastacmd -d "+sys.argv[1]+" -p F -s "+chr+" -L "+str(start_pos)+','+str(end_pos))
				print (str(pos)+' '+str(pos+len(seq)- 1)+' '+str(chr))
				len_rep = 0
				pos += (len(seq))
				pos += len(token[len(token)-1])
				seq = ''
				gap_found=False
			elif(line[len(line)-1] == 'N'):
				len_rep += len(line)
				seq += line
		elif(len(token)==1):
			pos += len(line)