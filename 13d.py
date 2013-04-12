### for removing the rep element and replacing it with N's

import sys

seq=''
for line in open(sys.argv[1],'r'):
	if(line[0]=='>'):
		print (line)
	else:
		line=line.strip()
		seq += line
		
		
nb=int(sys.argv[2])

rep=seq[nb:(len(seq)-nb)]

def replacedN(sequence):
	  complement = {'A':'N','C':'N','G':'N','T':'N','N':'N'}
	  return "".join([complement.get(nt.upper(), '') for nt in sequence[::-1]])

rep = replacedN(rep)



print (seq[0:nb]+rep+seq[(len(seq)-nb):])