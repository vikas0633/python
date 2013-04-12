#12n.py - for normalizing profiles

import sys
import numpy as np
sum_lib=[]
### load lib abundance
first_line=True
for line in open(sys.argv[1],'r'):
	line=line.strip()
	if first_line == False:
		if(len(line)>0):
			token=line.split('\t')
			for i in range(1,len(token)):
				sum_lib[i-1] += int(token[i]) 
				 
	else:
		token=line.split('\t')
		for i in range(0,len(token)):
			sum_lib.append(0)
		first_line = False
	
### take out lowest library count
minimum=min(sum_lib)

### normalization factor
sum_lib = [minimum/float(i) for i in sum_lib]

for line in open(sys.argv[1],'r'):
	line = line.strip()
	row=np.array([])
	if(len(line)>0):
		if(line[0]=='#'):
			header=line+'\t'+'sum'
			tokens=header.split('\t')
			for i in range(1,len(tokens)):
				header += ('\t'+tokens[i]+"_norm")
			header += "\t"+'score'
			print header
		else:
			tokens=line.split('\t')
			for i in range(1,len(tokens)-1): ### modify to have sum_norm as sum of all the other norm columns 
				if(float(tokens[i])==0):
					tokens[i]=1
				value=float(tokens[i])*float(sum_lib[i-1])
				line += '\t'+str(value)
				row=np.append(row,value)
			line += '\t'+str(sum(row))
			print(line+'\t'+str(np.std(row)/np.sqrt(np.average(row))))