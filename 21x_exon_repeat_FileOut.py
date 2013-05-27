#21x_exon_repeat.py- find the exon Repeat over lap - /Users/vikas0633/Desktop/script/python
# idea is to hash the repeat region and get the overlap with exonic regions
# input should be two gff3 files


# Usage: python ~/script/python/21x_exon_repeat.py -i RepeatMasker.gff -g sample.gff3 

import os,sys,getopt, re


### main argument to 

def options(argv):
    inputfile = ''
    gfffile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:g:o:",["ifile=","gfile=","out=",])
    except getopt.GetoptError:
       print 'python 21v_format_gff3.py -i <inputfile> -g <gfile> -o <out>'
       sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print 'python 21v_format_gff3.py -i <inputfile> -g <gfile> -o <out>'
    		sys.exit()
        elif opt in ("-i", "--ifile"):
        	inputfile = arg
    	elif opt in ("-g", "--gfile"):
          gfffile = arg
    	elif opt in ("-o", "--out"):
          out = arg
    
    return inputfile, gfffile, out

def get_size(file):
	size = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0:
			if line[0] != '#':
				token = line.split('\t')
				
				if token[0] not in size:
					size[token[0]] = int(token[4])
					
				else:
					if int(token[4]) > size[token[0]]:
						size[token[0]] = int(token[4])			
	return size

def hash_repeat(file,chr):
	hash = {}
	for line in open(file,'r'):
		line = line.strip()
		token = line.split('\t')
		
		if token[0] == chr:
			for i in range(int(token[3]),int(token[4])+1):
				hash[i] = ''
	return hash
		
def find_repetative_exons(file,hash,chr,out):

	hashSet = set(hash)
	
	first_line = True
	
	def printResult(ID,exon_count,exon_len,repeat_overlap):
        
		out.write(ID+'\t'+str(exon_count)+'\t'+str(exon_len)+'\t'+str(repeat_overlap)+'\t'+str(round(repeat_overlap/float(exon_len)*100,2))+'\n')
		
	ID = ''
	for line in open(file,'r'):
		line = line.strip()
		token = line.split('\t')
		if len(token) > 3:
			if token[0] == chr:
				if token[2]== 'mRNA':
					if first_line == False:
						printResult(ID,exon_count,exon_len,repeat_overlap)
					first_line = False
					exon_count = 0
					exon_len = 0
					repeat_overlap = 0
					### get the mRNA ID
					match = re.search(r'ID=.+;',line)
					if match:
						match = match.group().split(';')[0].replace('ID=','')
						ID = match	
				
				if token[2] == 'exon':
					exon_count += 1
					exon_len += int(token[4]) - int(token[3]) + 1
					
					### hash exonic positions
					temp = {}
					for i in range(int(token[3]),int(token[4])+1):
						temp[i] = ''
					tempSet = set(temp)
					
					### add the ovelap
					repeat_overlap += len(tempSet.intersection(hashSet))
	if ID != '':
		printResult(ID,exon_count,exon_len,repeat_overlap)
		
if __name__ == "__main__":
    
    file, geneModel, out = options(sys.argv[1:])
    
    ### open outfile
    outfile = open(out,'w')
    
    ### get the maximum chromosome size
    size = get_size(geneModel)
    
    ### print header
    outfile.write('ID'+'\t'+'Exon count'+'\t'+'Exon length'+'\t'+'Repeat overlap'+'\t'+'Repeat fraction'+'\n')
    
    ### hash the repeat file co-ordinates
    for chr in sorted(size.keys()):
    	hash = hash_repeat(file,chr)
    
    	### find exonic overlap for each mRNA
    	find_repetative_exons(geneModel,hash,chr,outfile)
    