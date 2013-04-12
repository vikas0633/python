### script to count the uniq MS supported genes

from B_hash_mRNA_IDs import *
import sys

###get the hash IDs
hash = hash_mRNA(sys.argv[1])


### open the protein file and count the gene hits

for line in open(sys.argv[2],'r'):
	line = line.strip()
	token = line.split('\t')
	
	if token[0] in hash:
		print token[0] 