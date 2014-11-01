#12n.py - for normalizing profiles

import sys
import numpy as np
sum_lib=[]
### load lib abundance
first_line=True
for line in open(sys.argv[1],'r'):
	line=line.strip()
	token=line.split('\t')
				 
	if first_line == True:
		for i in range(1,len(token)):
			sum_lib.append(0)
		first_line = False
	else:
		if(len(line)>0):
			for i in range(1,len(token)):
				sum_lib[i-1] += float(token[i]) 
	
### take out lowest library count
minimum=min(sum_lib)

## testing
'''
for lib in sum_lib:
	print str(sum_lib.index(lib)) +'\t'+ str(lib)
'''

### normalization factor
sum_lib = [minimum/float(i) for i in sum_lib]

first_line=True
for line in open(sys.argv[1],'r'):
	line = line.strip()
	row=np.array([])
	if(len(line)>0):
		if first_line == True:
			header=line+'\t'+'sum'
			tokens=header.split('\t')
			for i in range(1,len(tokens)):
				header += ('\t'+tokens[i]+"_norm")
			header += "\t"+'score'
			print header
			first_line = False
		else:
			tokens=line.split('\t')
			line += '\t'+str(sum([float(n) for n in tokens[1:]]))
			for i in range(1,len(tokens)): ### modify to have sum_norm as sum of all the other norm columns 
				if(float(tokens[i])==0):
					tokens[i]=1
				value=round(float(tokens[i])*float(sum_lib[i-1]),3)
				line += '\t'+str(value)
				row=np.append(row,value)
			line += '\t'+str(sum(row))
			print(line+'\t'+str(round(np.std(row)/np.sqrt(np.average(row)),3)))