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



def parse(file,chr,hash,key,print_gene_models,gene_id):
	for line in open(file,'r'):
		line = line.strip()
		token = line.split('\t')
		if len(line) > 0:		
			if line[0] != '#': 
				if (token[1] == key) & (token[0] == chr) & (token[2] == "gene"):
					for i in range(int(token[3]),int(token[4])+1):
						hash[i]=''
					if print_gene_models == True:
						if (token[1] == "CUFFLINKS") or  (token[1] == "GlimmerHMM") or (token[1] == "Augustus") or (token[1] == "GeneMark.hmm"):
							ID=(line.split('=')[1]).split(';')[0]
							gene_id[ID]=''
							print line
						elif (token[1] == "Ljchr0-6_pseudomol_20120830.chlo.mito.fa"):
							ID=(line.split('=')[1]).split(';')[0]
							gene_id[ID]=''
							mRNA_id = ID.replace('path','mrna')
							gene_id[mRNA_id]=''
							print line
						
						else:
							ID=(line.split('=')[2]).split(';')[0]
							gene_id[ID]=''
							print line
						
	return hash,gene_id

def filter_out(chr,hash,file):
	o = open('temp','w')
	for line in open(file,'r'):
		temp={}
		found = False
		line = line.strip()
		token = line.split('\t')
		if token[0] == chr:
			for i in range(int(token[3]),int(token[4])+1):
				if i in hash:
					found = True
					break
			if found == False:
				o.write(line+'\n') 
	o.close()


def find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models,gene_id):
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
				temp_set = set(temp)
				if (len(temp_set.intersection(hash_cov_set)) > 0.8*(len(temp_set)))  or (len(temp_set.intersection(hash_aug_set)) > 0.8*(len(temp_set))):
					if (token[1] == "CUFFLINKS") or  (token[1] == "GlimmerHMM") or (token[1] == "Augustus") or (token[1] == "GeneMark.hmm"):
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID]=''
						print line
					
					elif (token[1] == "Ljchr0-6_pseudomol_20120830.chlo.mito.fa"):
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID]=''
						mRNA_id = ID.replace('path','mrna')
						gene_id[mRNA_id]=''
						print line
						
					else:
						ID=(line.split('=')[2]).split(';')[0]
						gene_id[ID]=''
						print line
					
					for i in range(int(token[3]),int(token[4])+1):
						hash[i]=''
					
	return hash,gene_id
		
def find_overlap_denovo(file,chr,hash_cov,hash_TC,key,print_gene_models,gene_id):
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
				temp_set = set(temp)
				if (len(temp_set.intersection(hash_cov_set)) > 0.8*(len(temp_set)))  or (len(temp_set.intersection(hash_TC_set)) > 10):
					if (token[1] == "CUFFLINKS") or  (token[1] == "GlimmerHMM") or (token[1] == "Augustus") or (token[1] == "GeneMark.hmm"):
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID]=''
						print line
					elif (token[1] == "Ljchr0-6_pseudomol_20120830.chlo.mito.fa"):
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID]=''
						mRNA_id = ID.replace('path','mrna')
						gene_id[mRNA_id]=''
						print line
					else:
						ID=(line.split('=')[2]).split(';')[0]
						gene_id[ID]=''
						print line
					for i in range(int(token[3]),int(token[4])+1):
						hash[i]=''
	return hash,gene_id
	
def gene_stru(gene_id,infile):
	
	### process infile
	for line in open(infile,'r'):
		line = line.strip()
		token = line.split('\t')
		if len(line) > 0:
			if (line[0] != '#'):
				if (token[2] != 'gene'):
					if (token[1] == "CUFFLINKS"):
						keys = (line.split('=')[-1]).split('.')
						key = keys[0]+'.'+keys[1]
						if key in gene_id:
							print line
					if (token[1] == 'Ljchr0-6_pseudomol_20120830.chlo.mito.fa'):
						keys = ((line.split('=')[3])).split(';')[0].split('.')
						key = keys[0]+'.'+keys[1]+'.'+keys[2]
						if key in gene_id:
							print line
					if (token[1] == 'Lotus_TC'):
						keys = (line.split('=')[3]).split('.')
						key = keys[0]
						if key in gene_id:
							print line
					if (token[1] == "Augustus"):
						keys = (line.split('=')[-1]).split('.')
						key = 'gene.'+keys[1]+'.'+keys[2].strip()
						if key in gene_id:
							print line
					if (token[1] == "GlimmerHMM"):
						keys = (line.split('=')[-1]).split('.')
						key = 'gene.'+keys[1]+'.'+keys[2]+'.'+keys[3].strip()
						if key in gene_id:
							print line
					if (token[1] == "GeneMark.hmm"):
						key = (line.split('=')[-1]).replace('t','g')
						if key in gene_id:
							print line
				
			

def call_gff(file):
	infile = file
	for chr in size:
	
		### store gene_ids
		gene_id = {}
		
		
		### hash coverage file
		hash={}
		print_gene_models = False
		key = 'cufflinks'
		hash_cov,gene_id = parse(file,chr,hash,key,print_gene_models,gene_id)
		
		### hash Augustus
		hash={}
		print_gene_models = False
		key = 'Augustus'
		hash_aug,gene_id = parse(file,chr,hash,key,print_gene_models,gene_id)	
		

		### hash TC
		hash={}
		print_gene_models = False
		key = 'Lotus_TC'
		hash_TC,gene_id = parse(file,chr,hash,key,print_gene_models,gene_id)
		
		
		### hash the cufflinks gff co-ordinates
		hash={}
		key = 'CUFFLINKS'
		print_gene_models = True
		hash,gene_id = parse(file,chr,hash,key,print_gene_models,gene_id) 
		filter_out(chr,hash,file)		
		
		### hash the velvet co-ordinates and find overlap against coverage file and Augustus
		print_gene_models = True
		file = 'temp'
		key = 'Ljchr0-6_pseudomol_20120830.chlo.mito.fa'
		hash = {}
		hash,gene_id = find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models,gene_id)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		### de-novo
		### hash the Augustus co-ordinates and find overlap against coverage file and TC
		print_gene_models = True
		file = 'temp'
		key = 'Augustus'
		hash = {}
		hash,gene_id = find_overlap_denovo(file,chr,hash_cov,hash_TC,key,print_gene_models,gene_id)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		
		### hash the GeneMark co-ordinates and find overlap against coverage file and TC
		print_gene_models = True
		file = 'temp'
		key = 'GeneMark.hmm'
		hash = {}
		hash,gene_id = find_overlap_denovo(file,chr,hash_cov,hash_TC,key,print_gene_models,gene_id)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		
		### hash the GlimmerHMM co-ordinates and find overlap against coverage file and TC
		print_gene_models = True
		file = 'temp'
		key = 'GlimmerHMM'
		hash = {}
		hash,gene_id = find_overlap_denovo(file,chr,hash_cov,hash_TC,key,print_gene_models,gene_id)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		
		### hash the TC co-ordinates and find overlap against coverage file and Augustus
		print_gene_models = True
		file = 'temp'
		key = 'Lotus_TC'
		hash = {}
		hash,gene_id = find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models,gene_id)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		'''
		### Homology based
		print_gene_models = True
		file = 'temp'
		key = 'Lotus_TC'
		hash = {}
		hash,gene_id = find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models,gene_id)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		'''
		
		
		
	
		### print gene structure
		gene_stru(gene_id,infile)
		
		file = infile

		
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
	
	
	