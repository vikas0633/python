### script for making igv files for clustering and visualization

import sys

def load_expression(infile):
	### open profile file and make a hash for storing information into the hash
	expression={}
	for line in open(infile,'r'):
		if(len(line)>0):
			line=line.strip()
			token=line.split('\t')
			expression[token[0]]=token[len(token)/2]
			for i in range(len(token)/2+1,len(token)):
				expression[token[0]] +='\t'+token[i]
	return expression
	
	
	
### open fasta file and 
def load_fasta(infile):
	fasta={}
	for line in open(infile,'r'):
		if(len(line)>0):
			line=line.strip()
			if(line[0]=='>'):
				key=line[1:]
			else:
				fasta[key]=line
	return fasta


### make igv files
def make_igv(sam_file,expression,fasta):
	print 'chromosome\tstart\tstop\tsequence\t'+expression['#Sequence']
	for line in open(sam_file,'r'):
		line=line.strip()
		if(len(line)>0):
			token=line.split('\t')
			id=token[0]
			chro=token[2]
			map=token[3]
			print chro+'\t'+map+'\t'+str(int(map)+int(len(fasta[id])))+'\t'+(token[0].split('_'))[1]+'\t'+expression[fasta[id]]

if __name__ == "__main__":
	
	### get the hash for expression values
	expression=load_expression(sys.argv[1])
	
	### get the hash for fasta sequences
	fasta=load_fasta(sys.argv[2])
	
	### run loop for sam file
	make_igv(sys.argv[3],expression,fasta)