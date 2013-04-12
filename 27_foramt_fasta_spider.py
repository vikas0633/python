
### 27_foramt_fasta_spider.py - /Users/vikas0633/Desktop/script/python/ - script to format the fasta headers according to the Thomas's explanations


### requirements:
'''
Here is the fasta header structure we discussed.
 
The header structure from the original files look like this:
>L9107_T1/1_Tarantula_S,L9107_T1/1_Tarantula_S
>L9109_T1/1_Tarantula_S,gi|307210400|gb|EFN86965.1|,1e-12,76.6, VWFA and cache domain-containing protein 1 [Harpegnathos saltator]
 
We would like to have the new structure look like this:
>L9107_T1/1_Tarantula_S_fr1,L9107_T1/1_Tarantula_S
>L9109_T1/1_Tarantula_S_fr3,VWFA and cache domain-containing protein 1 [Harpegnathos saltator], gi|307210400|gb|EFN86965.1|,1e-12,76.6
 
The changes are that _frX is added to the accession and that the name of the
protein is moved to the front (the name should be everything after the 4th
comma in the original header). As mentioned then we use everything before the
first comma as the accession for the protein and everything after as the
protein name. 
'''

### check output with 
# python ~/script/python/27_foramt_fasta_spider.py -i  

import os,sys,getopt, re
from C_loadFasta import *

### main argument to 

def options(argv):
	inputfile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:",["ifile="])
	except getopt.GetoptError:
		print 'python 27_foramt_fasta_spider.py -i <inputfile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'python 27_foramt_fasta_spider.py -i <inputfile>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
			
	return inputfile
    
							
def format_fasta(fasta):
	for line in open(fasta,'r'):
		line = line.strip()
		if len(line) > 0:
			if line[0] != '#':
				if line[0] == '>':
					if re.search('UNIQ',line):
						token = line[1:].split('UNIQ') ### remove the annotations aded by the translateDNA.py script
					elif re.search(' DUPL_LONG ',line):
						token = line[1:].split('DUPL_LONG')
					### case where no blast annotation is present
					tokens = token[0].split(',')
					if len(tokens)==2: 							#13_T2/11_Tarantula_WB,L13_T2/11_Tarantula_WB_fr4 
						print '>'+tokens[1].strip()+','+tokens[0]
					else:										# >L14_T4/4_Tarantula_WB,gi|16930529|gb|AAL31950.1|,4e-40,_159, CDH1-D [Gallus gallus]_fr6
						print '>' + tokens[0]+token[0][-5:-1] + ',' + tokens[-1][0:len(tokens[-1])-4].replace('_','').lstrip()+',' + (','.join(tokens[1:4])).replace('_','')
				else:
					print line
		

if __name__ == "__main__":
	fasta = options(sys.argv[1:])

	format_fasta(fasta)


