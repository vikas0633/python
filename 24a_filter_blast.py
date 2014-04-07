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
		query = token[0]

		annotation = hash[token[1]][1:] 
		identity = token[2]
		e_value = token[10]
		key1 = 'Uncharacterized'
		key2 = 'Putative'
		key3 = 'Hypothetical'
		key4 = 'novel'
		key5 = 'Unknown'
		
		
		
		if first_line == True:
			x = 1
		if query == last_query:
			if  re.search(key1.lower(),last_annotation.lower()):
				continue
			elif re.search(key2.lower(),last_annotation.lower()):
				continue
			elif re.search(key3.lower(),last_annotation.lower()):
				continue
			elif re.search(key4.lower(),last_annotation.lower()):
				continue
			elif re.search(key5.lower(),last_annotation.lower()):
				continue
			else:
				if flag == False:
					print last_query +'\t'+ last_annotation +'\t'+ last_identity +'\t'+ last_e_value
					flag = True
		
		if (query != last_query) & (first_line == False): ### check for new query
			
			if flag == False:
				if re.search(key1.lower(),last_annotation.lower()):
					print last_query +'\t'+'Non Chatacterized Hit- '+last_annotation +'\t'+ last_identity +'\t'+ last_e_value ### print for last query if all were putative
				elif re.search(key2.lower(),last_annotation.lower()):
					print last_query +'\t'+'Non Chatacterized Hit- '+last_annotation +'\t'+ last_identity +'\t'+ last_e_value ### print for last query if all were putative
				else:
					print last_query +'\t'+ last_annotation +'\t'+ last_identity +'\t'+ last_e_value ### print for last query if all were putative
			flag = False
		
		last_query = query
		last_annotation = annotation 
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
