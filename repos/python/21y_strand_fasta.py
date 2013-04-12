## 21y_strand_fasta.py - script takes a GFF3 file and correct fasta file if minus strand - /Users/vikas0633/Desktop/script/python

## this script was tested using Augustus GFF3 file

import os,sys,getopt, re

### Usage
'''
python ~/script/python/21y_strand_fasta.py -f CDS_20121026_only_exon_augustus.fa -g 20121026_only_exon_augustus.gff3
'''


### main argument to 

def options(argv):
    file = ''
    identifier = ''
    try:
       opts, args = getopt.getopt(argv,"hf:g:",["fasta=","gff3="])
    except getopt.GetoptError:
    	print '''
    				python 21y_strand_fasta.py 
    					-f <fasta>
    					-g <gff3> 
    			'''       
    	sys.exit(2)
    
    for opt, arg in opts:
    	if opt == '-h':
			print '''
						python 21y_strand_fasta.py 
							-f <fasta>
							-g <gff3> 
					'''       
			sys.exit()
        elif opt in ("-f", "--fasta"):
        	fasta = arg
    	elif opt in ("-g", "--gff3"):
        	gff3 = arg
        	
    return fasta, gff3
    
### function for reverse complement
def reverseComplement(sequence):
  complement = {'A':'T','C':'G','G':'C','T':'A','N':'N'}
  return "".join([complement.get(nt.upper(), '') for nt in sequence[::-1]])


	
'''
scaffold5730	AUGUSTUS	CDS	1	226	0.56	-	0	ID=g1.t1.cds;Parent=g1.t1
scaffold6210	AUGUSTUS	CDS	2540	2818	0.86	+	0	ID=g2.t1.cds;Parent=g2.t1
'''

def hash_gff3(gff3):
	hash = {}
	''' This function take a gff3 file and stores ID as key and strand as value'''
	for line in open(gff3,'r'):
		line = line.strip()
		token = line.split('\t') 	
		match = re.search(r'ID=.+;',line)
		ID = match.group().split(';')[0].replace('ID=','')
		hash[ID] = token[6]
	
	return hash
	
### reverse complement the negative strand
def process(fasta,hash):
	first_line = True
	for line in open(fasta,'r'):
		line = line.strip()
		if line[0]=='>':
			if first_line == False:
				print '>'+header
				### test for strand
				if hash[header]=='-':
					print reverseComplement(seq)
				else:
					print seq
			
			header = line[1:].split(' ')[0]
			seq = ''
		else:
			seq += line
		
		first_line = False
	
	### for last fasta file
	print '>'+header
	### test for strand
	if hash[header]=='-':
		print reverseComplement(seq)
	else:
		print seq
	
    
if __name__ == "__main__":
    
    fasta, gff3 = options(sys.argv[1:])
    
    ### fetch strands from the gff3 file
    hash = hash_gff3(gff3)
    
    ### process the fasta file
    process(fasta,hash)
    
    