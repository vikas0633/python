#21v_format_gff3.py- script to format gff3 file in order to put in MySQL table - /Users/vikas0633/Desktop/script/python

# Usage: python /Users/vikas0633/Desktop/script/python/21v_format_gff3.py -i sample.gff3

import os,sys,getopt, re


### main argument to 

def options(argv):
    inputfile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:",["ifile=",])
    except getopt.GetoptError:
       print 'python 21v_format_gff3.py -i <inputfile> '
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print 'python 21v_format_gff3.py -i <inputfile> '
          sys.exit()
       elif opt in ("-i", "--ifile"):
          inputfile = arg
    
    return inputfile

def load_elements(inf):
	elements = {}
	for line in open(inf,'r'):
		line = line.strip()
		### find parent name
		token = line.split('\t')
		if len(token) > 3:
			if token[2] == "mRNA":
				if (token[1]=="CUFFLINKS" or token[1]=="Cufflinks"):
					match = re.search(r'ID=.+;',line)
					if match:
						match = match.group().split(';')[0].replace('ID=','')
						mRNA = match
						if match in elements:
							elements[match] += line+'\n'
						else:
							elements[match] = line+'\n'
				else:
					match = re.search(r'Parent=.*',line)
					if match:
						match = match.group().split(';')[0].replace('Parent=','')
						mRNA = match
						if match in elements:
							elements[match] += line+'\n'
						else:
							elements[match] = line+'\n'
			if (token[2] != "gene") & (token[2] != "mRNA"):
				elements[mRNA] += line+'\n'
	return elements
	
def write_elements(inf,elements):
	for line in open(inf,'r'):
		line = line.strip()
		### find parent name
		token = line.split('\t')
		if len(token) > 3:
			if token[2] == "gene":
				if (token[1]=="CUFFLINKS" or token[1]=="Cufflinks"):
					match = re.search(r'Name=.+',line)
					match = match.group().split(';')[0].replace('Name=','')
					if match in elements:
						print line
						print elements[match]
				else:
					match = re.search(r'ID=.+;',line)
					if match:
						match = match.group().split(';')[0].replace(';','').replace('ID=','') 
						if match in elements:
							print line
							print elements[match]

if __name__ == "__main__":
    
    inf = options(sys.argv[1:])

    ### parse file to get the gene names
    elements = load_elements(inf)
    
    ### use the gene names to find other genomic elements
    write_elements(inf,elements)
