### 21q_combine_GTF.py
### This is the script for combining various annotations files

### all evidences must be unique

### python ~/Desktop/script/python/21q_combine_GTF.py sample.05_transcripts.MG20_mRNA_illumina_denovo.gff3

### sample gff
'''
Ljchr1_pseudomol_20120830	GlimmerHMM	mRNA	2404	3747	.	-	.	ID=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudomol_20120830.path1.gene1
Ljchr1_pseudomol_20120830	GlimmerHMM	CDS	2404	2514	.	-	0	ID=Ljchr1_pseudomol_20120830.cds1.1;Parent=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudomol_20120830.path1.gene1;Note=final-exon
Ljchr1_pseudomol_20120830	GlimmerHMM	CDS	2953	3108	.	-	0	ID=Ljchr1_pseudomol_20120830.cds1.2;Parent=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudomol_20120830.path1.gene1;Note=internal-exon
Ljchr1_pseudomol_20120830	GlimmerHMM	CDS	3198	3278	.	-	0	ID=Ljchr1_pseudomol_20120830.cds1.3;Parent=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudom
'''

import sys,os

import datetime
now = datetime.datetime.now()

#o = open(str(now.strftime("%Y-%m-%d_%H%M_"))+'gene.gff','w')

def options(argv):
	
	gff = argv
	
	return gff



def parse(file,chr,hash,key,print_gene_models):
	for line in open(file,'r'):
		line = line.strip()
		token = line.split('\t')
		if len(line) > 0:		
			if line[0] != '#': 
				if (token[1] == key) & (token[0] == chr) & (token[2] == "gene"):
					for i in range(int(token[3]),int(token[4])+1):
						hash[i]=''
					if print_gene_models == True:
						print line	
	return hash

def filter_out(chr,hash,file):
	o = open('temp','w')
	for line in open(file,'r'):
		line = line.strip()
		token = line.split('\t')
		if token[0] == chr:
			if int(token[3]) in hash:
				continue
			elif int(token[4]) in hash:
				continue
			elif (int(token[3])+int(token[4]))/2 in hash:
				continue
			else:
				o.write(line+'\n') 
	o.close()


def find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models):
	hash = {}
	hash_cov_set = set(hash_cov)
	hash_aug_set = set(hash_aug)
	for line in open(file,'r'):
		temp = {}
		line = line.strip()
		token = line.split('\t')
		if line[0] != '#': 
			if (token[1] == key) & (token[0] == chr) & (token[2] == "gene"):
				for i in range(int(token[3]),int(token[4])+1):
					temp[i]=''
					hash[i]=''
				temp_set = set(temp)
				if (len(temp_set.intersection(hash_cov_set)) > 0.8*(len(temp_set)))  or (len(temp_set.intersection(hash_aug_set)) > 0.8*(len(temp_set))):
					print line
					
	return hash
		
def find_overlap_denovo(file,chr,hash_cov,hash_TC,key,print_gene_models):
	hash = {}
	hash_cov_set = set(hash_cov)
	hash_TC_set = set(hash_TC)
	for line in open(file,'r'):
		temp = {}
		line = line.strip()
		token = line.split('\t')
		if line[0] != '#': 
			if (token[1] == key) & (token[0] == chr) & (token[2] == "gene"):
				for i in range(int(token[3]),int(token[4])+1):
					temp[i]=''
					hash[i]=''
				temp_set = set(temp)
				if (len(temp_set.intersection(hash_cov_set)) > 0.8*(len(temp_set)))  or (len(temp_set.intersection(hash_TC_set)) > 10):
					print line
					
	return hash

def call_gff(file):
	for chr in size:
		### hash the cufflinks gff co-ordinates
		hash={}
		key = 'CUFFLINKS'
		print_gene_models = True
		hash = parse(file,chr,hash,key,print_gene_models) 
		filter_out(chr,hash,file)
		
		### hash coverage file
		hash={}
		print_gene_models = False
		key = 'cufflinks'
		hash_cov = parse(file,chr,hash,key,print_gene_models)
		
		### hash Augustus
		hash={}
		print_gene_models = False
		key = 'Augustus'
		hash_aug = parse(file,chr,hash,key,print_gene_models)	
		

		### hash TC
		hash={}
		print_gene_models = False
		key = 'Lotus_TC'
		hash_TC = parse(file,chr,hash,key,print_gene_models)
		
		
		### hash the velvet co-ordinates and find overlap against coverage file and Augustus
		print_gene_models = True
		file = 'temp'
		key = 'Ljchr0-6_pseudomol_20120830.chlo.mito.fa'
		hash = {}
		hash = find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		### de-novo
		### hash the Augustus co-ordinates and find overlap against coverage file and TC
		print_gene_models = True
		file = 'temp'
		key = 'Augustus'
		hash = {}
		hash = find_overlap_denovo(file,chr,hash_cov,hash_TC,key,print_gene_models)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		### hash the GeneMark co-ordinates and find overlap against coverage file and TC
		print_gene_models = True
		file = 'temp'
		key = 'GeneMark.hmm'
		hash = {}
		hash = find_overlap(file,chr,hash_cov,hash_TC,key,print_gene_models)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		### hash the GlimmerHMM co-ordinates and find overlap against coverage file and TC
		print_gene_models = True
		file = 'temp'
		key = 'GlimmerHMM'
		hash = {}
		hash = find_overlap(file,chr,hash_cov,hash_TC,key,print_gene_models)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		
		### hash the TC co-ordinates and find overlap against coverage file and Augustus
		print_gene_models = True
		file = 'temp'
		key = 'Lotus_TC'
		hash = {}
		hash = find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
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

if __name__ == "__main__":
	
	file = options(sys.argv[1:])[0] 
	
	### get the end point for each chromosome
	size = get_size(file)
	
	### call gff file one by one
	call_gff(file)
	
	
	