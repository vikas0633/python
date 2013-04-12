#21ac_addType.py - script to add gene type - /Users/vikas0633/Desktop/script/python

# Usage: python ~/script/python/21ac_addType.py -i sample.gff3 -t type

import os,sys,getopt, re


### main argument to 

def options(argv):
	inputfile = ''
	type = ''
	try:
		opts, args = getopt.getopt(argv,"hi:t:",["ifile=","type="])
	except getopt.GetoptError:
		print 'python 21v_format_gff3.py -i <inputfile> -t <type>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'python 21v_format_gff3.py -i <inputfile> -t <type>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-t", "--type"):
			type = arg

	return inputfile, type


### hash the type
def HASH(type):
	hash = {}
	for line in open(type,'r'):
		if line[0]!='#':
			if len(line) > 2:
				line = line.strip()
				token = line.split('\t')
				if (token[0][0:4] != "CUFF") & (token[0][0:2]!='gi'): ### cufflinks gene models have .1 .2 for each gene while other methods only use .1, .2 for various mRNAs
					key = ".".join(token[0].split(".")[0:len(token[0].split("."))-1])
				else:
					key = token[0]
				key = key.replace('mrna','path') ### replace the gmap models so that match the gene name
				if len(token)==1:
					hash[key] = ' '
				else:
					hash[key] = token[1]
	return hash

### load mRNAids and respective models
def load_elements(inf,hash):
	elements = {}
	for line in open(inf,'r'):
		line = line.strip()
		### find parent name
		token = line.split('\t')
		if len(token) > 3:
			if line[-1]==';':
				line = line[0:len(line)-1]
			if token[2] == "mRNA":
				if token[1]=="CUFFLINKS":
					match = re.search(r'ID=.+;',line)
					if match:
						match = match.group().split(';')[0].replace('ID=','')
						mRNA = match
						match = match.replace('mrna','path')
						if match in elements:
							elements[match] += line + '; Type="'+hash[match]+'"\n'
						else:
							elements[match] = line + '; Type="'+hash[match]+'"\n'
				else:
					match = re.search(r'Parent=.*',line)
					if match:
						match = match.group().split(';')[0].replace('Parent=','')
						mRNA = match
						match = match.replace('mrna','path')
						
						### add the protein coding type for conserved proteins
						if match.startswith('gi'):
							hash[match] = 'protein_coding'
						
						if match in elements:
							elements[match] += line +'; Type="'+hash[match]+'"\n'
						else:
							elements[match] = line +'; Type="'+hash[match]+'"\n'
			if (token[2] != "gene") & (token[2] != "mRNA"):
				elements[mRNA] += line +'; Type="'+hash[match]+'"\n'
	return elements
	
def write_elements(inf,elements):
	for line in open(inf,'r'):
		line = line.strip()
		### find parent name
		token = line.split('\t')
		if len(token) > 3:
			if token[2] == "gene":
				if token[1]=="CUFFLINKS":
					match = line.split('%')[-1].replace(';','')[2:]
					if match in elements:
						print line + '; Type="'+hash[match]+'"'
						print elements[match]
				else:
					match = re.search(r'ID=.+;',line)
					if match:
						match = match.group().replace(';','').replace('ID=','') 
						if match in elements:
							print line + '; Type="'+hash[match]+'"'
							print elements[match]

if __name__ == "__main__":
    
	inf, type = options(sys.argv[1:])
	
	### hash the type values
	hash = HASH(type)


	### parse file to get the gene names
	elements = load_elements(inf,hash)
	
	### use the gene names to find other genomic elements
	write_elements(inf,elements)
