### 24a_filter_blast.py - /Users/vikas0633/Desktop/script/python/ - script for filtering blast results 

### input is two files, one with blast headers and other with blast hits

## Remove the following keywords:

## Putative
## Uncharacterized

# Usage: 
'''
python ~/Desktop/script/python/24a_filter_blast.py \
-i uniprot_sprot.fasta.headers \
-b Ljr_cds_protein.Ljr3.0.20130102.refined.fa.blastout.sorted \
> Ljr_cds_protein.Ljr3.0.20130102.refined.fa.blastout.genewise 

'''

import os,sys,getopt, re

### main argument to 

def options(argv):
    inputfile = ''
    blast = ''
    try:
    	opts, args = getopt.getopt(argv,"hi:b:",["headers=","blast="])
    except getopt.GetoptError:
    	print 'python 21v_format_gff3.py -i <headers> -b <blast>'
    	sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print 'python 21v_format_gff3.py -i <headers> -b <blast>'
    		sys.exit()
    	elif opt in ("-i", "--headers"):
    		inputfile = arg
    	elif opt in ("-b", "--blast"):
    		blast = arg
    
    return inputfile, blast
    
def hash_headers(file):
	headers = {}
	for line in open(file,'r'):
		line = line.strip()
		token = line.split(' ')
		headers[token[0][1:]] = line
	return headers
	
def filter_blast(hash,blast):
	last_query = ''
	first_line = True
	flag = False
	for line in open(blast,'r'):
	    line = line.strip()
	    token = line.split('\t')
	    query = token[0].split(',')[0]
	    flag2 = True 
	    annotation = hash[token[1]][1:].split(']')[0]+']'
	    annotation2 = ' '.join(annotation.split(' ')[1:]) + ' '+annotation.split(' ')[0]
	    identity = token[2]
	    e_value = token[10]
	
	    if first_line == False:
		if last_query != query:
		    if flag == False:
		        print last_query +'\t'+last_annotation +'\t'+ last_identity +'\t'+ last_e_value ### print for last query if all were putative
		    flag = False
	    
	    filter_key = ['Uncharacterized', 'Putative','Unnamed','Predicted','hypothetical protein', 'Uncharacterized protein', 'unknown']
	    for key in filter_key:
		if re.search(key.lower(), annotation2):
		    flag2 = False 
		    break
	    if flag == False:
		if flag2 == True:
		    print query +'\t'+annotation2 +'\t'+ identity +'\t'+ e_value ### print for last query if all were putative
		    flag = True
		

	    
	    last_query = query
	    last_annotation = annotation2 
	    last_identity = identity
	    last_e_value = e_value 
	    first_line = False
		
	if flag == False:
	    print last_query +'\t'+last_annotation +'\t'+ last_identity +'\t'+ last_e_value ### print for last query if all were putative
		
		
if __name__ == "__main__":
    
    headers, blast = options(sys.argv[1:])
    
    ### hash the headers using space as divider
    hash = hash_headers(headers)
    
    ### filter the blast results
    filter_blast(hash,blast)
