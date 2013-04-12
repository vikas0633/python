###21d_take_out_gene.py - this script takes out a sequence from fasta file given correct header name 
### nice -n 19 python 21d_take_out_gene.py -f <fasta file> -n <header_name>
### python 21d_take_out_gene.py -f Kasuza.fasta.refined -n chr5_31999127_31999213

import sys, getopt, os

def file_empty(file):
    count = sum([1 for line in open(file)])
    if count == 0:
        sys.exit(file+' is empty')

### get the options 
def options(argv):
	complete = False; GTF = False; start = ''; end = '';file=''
	try:
		opts, args = getopt.getopt(argv,"hf:n:s:e:g",["fasta=","header=","start=","end=","GTF="])
	except getopt.GetoptError:
		print	''' 
			Usages:
				python 21d_take_out_gene.py 
						-f <fasta file> 
						-n <header_name> 
						-s <start>  [leave empty if whole sequence need to be extracted] 
						-e <end>
						-g <GTF> [by defualt fasta file]
				'''
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print	''' 
			Usages:
				python 21d_take_out_gene.py 
						-f <fasta file> 
						-n <header_name> 
						-s <start>  [leave empty if whole sequence need to be extracted] 
						-e <end>
						-g <GTF> [by defualt fasta file]
				'''
			sys.exit()
		elif opt in ("-f", "--fasta"):
			file = arg
		elif opt in ("-n", "--header"):
			header = arg
		elif opt in ("-s", "--start"):
			start = int(arg)
		elif opt in ("-e", "--end"):
			end = int(arg)
		elif opt in ("-g", "--GTF"):
			GTF = True		
	return file, header, start, end , GTF, complete

def extract_seq(file, header, start, end , GTF, complete):
	header = header.replace(' ','').replace('>','')
	flag = False
	length = 0
	seq = ''
	if complete == False:
		os.system('nice -n 19 fastacmd -d '+file+' -p F -s '+header+' -L '+str(start)+','+str(end))
	else:
		for line in open(file,'r'):
			if len(line) > 0:
				if line[0] != '#':
					if line[0]=='>':
						line = line.strip()
						if line[1:] == header:
							flag = True
							print line
						else:
							flag = False
					else:
						if flag == True:
							line = line.strip()
							print line

def extract_gtf(file, header, start, end , GTF, complete):
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0:
			if line[0] != '#':
				token = line.split('\t')
				if complete == True:
					if token[0] == header:
						print line
				else:
					if (token[0] == header) & (start <= int(token[3])) & (int(token[4]) <= end):
						print line 
	
if __name__ == "__main__":
    
    # get the file and header
    file, header, start, end , GTF, complete = options(sys.argv[1:])
    if start == '':
       		complete = True
    
    file_empty(file)
    
    index = True
    
    if index == False:
    	os.system('formatdb -i '+file+' -p F -o T')
    if GTF == True:
    	extract_gtf(file, header, start, end , GTF, complete)
    else:
		# get the sequence
		extract_seq(file, header, start, end , GTF, complete)
			
    
