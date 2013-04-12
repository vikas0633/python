### 25a_plot_gene_freq.py - /Users/vikas0633/Desktop/script/python/ - script for plotting gene frequencies across each chromosome

### idea is to find the place with the high gene density and write the chromosomal positions

### take a gff file containing gene locations
### 3rd column must have 'gene' keyword

# Usage: python ~/script/python/25a_plot_gene_freq.py -i 20121227_conserved_Ljr3_0.refined.gff3

import os,sys,getopt, re
import numpy as np


### main argument to 

def options(argv):
    inputfile = ''
    window = 1000
    step = 100
    
    try:
       opts, args = getopt.getopt(argv,"hi:w:s:",["gfffile=","window=","step="])
    except getopt.GetoptError:
       print 'python 25a_plot_gene_freq.py -i <gfffile> -w <window> -s <step>'
       sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print 'python 25a_plot_gene_freq.py -i <gfffile> -w <window> -s <step>'
    		sys.exit()
        elif opt in ("-i", "--ifile"):
        	infile = arg
    	elif opt in ("-w", "--window"):
        	window = arg
        elif opt in ("-s", "--step"):
        	step = arg
    	
    return infile, window, step

def get_chro(file):
	'''this function return a unique list of chromosome/contigs name '''
	hash = {}
	for line in open(file,'r'):
		line = line.strip()
		token = line.split('\t')
		if len(token) > 3:
			if token[0] in hash:
				if int(token[3]) > hash[token[0]]:
					hash[token[0]] = int(token[3])
			else:
				hash[token[0]] = int(token[3])
	return hash
    
    
def find_gene_density(file, window, step):
	''' function to find gene density in a given window '''
	chro = get_chro(file)
	for key in sorted(chro):
		list = []
		gene = {}
		for line in open(file,'r'):
			line = line.strip()
			token = line.split('\t')
			if len(token) > 3:
				if (token[2] == "gene") & (token[0] == key):
					### note gene position in a list
					if int(token[3]) in gene:
						gene[int(token[3])] += 1
					else:
						gene[int(token[3])] = 1
		''' 
		out of memory error
		### hash the co-ordinates
		pos = {}
		for i in range(1,chro[key]+1):
			if i in gene:
				pos[i] = gene[i]
			else:
				pos[i] = 0
		'''
		
		### check the average start for each chromosome
		#print sum(pos.values())/float(len(pos.values()))
		
		### let's scan through the genome
		last_start = 0
		gene_count = 0
		temp = 0
		flag = True
		for i in range(window/2,chro[key]-window/2,step):
			genes = 0
			for j in range(i-window/2,i+window/2,1):
				if j in gene:
					genes += gene[j]
			if genes > 5:
				if (i - last_start > step) & (gene_count > 0):
					print key, start, last_start, gene_count
					temp = 0	
					gene_count = 0
					flag = True
				if flag == True:
					start = i
					flag = False			
				last_start = i
				gene_count += abs(temp - genes) 
				temp = genes
				
		## plot_gene_density
		
		#list = np.array(list)
		#print np.average(list)
	
    
if __name__ == "__main__":
    
    file, window, step = options(sys.argv[1:])
    
    find_gene_density(file, window, step)
    
