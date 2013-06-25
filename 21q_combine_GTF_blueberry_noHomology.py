### 21q_combine_GTF.py
### This is the script for combining various annotations files

### all evidences must be unique

### python ~/Desktop/script/python/21q_combine_GTF.py sample.20121108_combined.gff3

### sample gff
'''
Ljchr1_pseudomol_20120830	GlimmerHMM	mRNA	2404	3747	.	-	.	ID=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudomol_20120830.path1.gene1
Ljchr1_pseudomol_20120830	GlimmerHMM	CDS	2404	2514	.	-	0	ID=Ljchr1_pseudomol_20120830.cds1.1;Parent=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudomol_20120830.path1.gene1;Note=final-exon
Ljchr1_pseudomol_20120830	GlimmerHMM	CDS	2953	3108	.	-	0	ID=Ljchr1_pseudomol_20120830.cds1.2;Parent=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudomol_20120830.path1.gene1;Note=internal-exon
Ljchr1_pseudomol_20120830	GlimmerHMM	CDS	3198	3278	.	-	0	ID=Ljchr1_pseudomol_20120830.cds1.3;Parent=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudom
'''

import sys,os, commands
import datetime
now = datetime.datetime.now()


### Global variable
#key is the second column in the GFF3/GTF file
COVERAGE = "cufflinks"  ### key for the file that contains the RNA-seq pileup in GTF
AUGUSTUS = "Augustus"
GlimmerHMM = "GlimmerHMM"
CUFFLINKS = "CUFFLINKS"
GENEMARK = "GeneMark.hmm"
Key_454 = "454Scaffolds"

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
				if (token[1] == key) & (token[0] == chr) & ((token[2] == "gene") or (token[2] == "exon")):
					for i in range(int(token[3]),int(token[4])+1):
						hash[i]=''
					if print_gene_models == True:
						if (token[1] == CUFFLINKS) or  (token[1] == GlimmerHMM) or (token[1] == AUGUSTUS) or (token[1] == GENEMARK):
							ID=(line.split('=')[1]).split(';')[0]
							gene_id[ID]=''
							print line
						elif (token[1] == Key_454):
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
					if (token[1] == CUFFLINKS) or  (token[1] == GlimmerHMM) or (token[1] == AUGUSTUS) or (token[1] == GENEMARK):
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID]=''
						print line
					
					elif (token[1] == Key_454):
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
					if (token[1] == CUFFLINKS) or  (token[1] == GlimmerHMM) or (token[1] == AUGUSTUS) or (token[1] == GENEMARK):
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID]=''
						print line
					elif (token[1] == Key_454):
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

def do_blast(file,chr,hash,print_gene_models,gene_id,query):
	flag = False
	### do the blast
	os.system('nice -n 19 blastx -db '+database+' -query temp.fa -evalue 0.0000000001 -outfmt 6 -out temp.blastout')
	### get the length of the query
	len_query=int((commands.getoutput('wc -c temp.fa')).split()[0]) - int((commands.getoutput('wc -l temp.fa')).split()[0]) - 5
	count = sum([1 for line in open('temp.blastout')])
	
	if count !=0:
		for line in open('temp.blastout','r'):
			line = line.strip()
			token = line.split()
			### check if query is aligned at-least 80%
			if abs(int(token[6]) - int(token[7])) > 0.8*(len_query): 
				flag = True
				break
		
	if flag == True:
		line = query.strip()
		token=line.split('\t')
		ID=(line.split('=')[1]).split(';')[0]
		gene_id[ID]=''
		print line+'.homology'
		for i in range(int(token[3]),int(token[4])+1):
			hash[i]=''
	return hash,gene_id
		
def homology_evidence(file,chr,hash,print_gene_models,gene_id,key):
	
	hash = {}
	### make blast database
	#formatdb -i Ljchr0-6_pseudomol_20120830.chlo.mito.fa -p F -o T
	#makeblastdb -in /Users/vikas0633/Desktop/temp/arabi_Brassi_lotus_medicago_protein.fasta
	new_gene=False
	os.system('echo ">"'+str('temp')+' >temp.fa')
	for line in open(file,'r'):
		flag = False
		line = line.strip()
		token = line.split('\t')
		if  token[0] == chr:
			if (token[1] == key):
				if (token[2]=="gene") :
					if new_gene == True:
						### get the fasta file
						if int((commands.getoutput('wc -c temp.fa')).split()[0]) > 10 : ### make sure that sequence is not empty [pileup gene models have no exons]
							hash,gene_id=do_blast(file,chr,hash,print_gene_models,gene_id,query)
						os.system('echo ">"'+str('temp')+' >temp.fa')
					new_gene = True
					query = line
				if token[2]=="exon":
					os.system('nice -n 19 fastacmd -d '+ref_seq+' -p F -s '+token[0]+' -L '+str(token[3])+','+str(token[4]) +' | awk "NR>1" >>temp.fa')
			
	
	return hash,gene_id
	
	

def gene_stru(gene_id,infile):
	
	### process infile
	for line in open(infile,'r'):
		line = line.strip()
		token = line.split('\t')
		if len(line) > 0:
			if (line[0] != '#'):
				if (token[2] != 'gene'):
					if (token[1] == CUFFLINKS):
						keys = (line.split('=')[-1]).split('.')
						key = keys[0]+'.'+keys[1]
						if key in gene_id:
							print line
					if (token[1] == Key_454):
						key = (line.split('=')[1]).split(';')[0]
						if key in gene_id:
							print line
					if (token[1] == AUGUSTUS):
						keys = (line.split('=')[-1]).split('.')
						key = 'gene.'+keys[1]+'.'+keys[2].strip()
						if key in gene_id:
							print line
					if (token[1] == GlimmerHMM):
						keys = (line.split('=')[-1]).split('.')
						key = 'gene.'+keys[1]+'.'+keys[2]+'.'+keys[3].strip()
						if key in gene_id:
							print line
					if (token[1] == GENEMARK):
						key = (line.split('=')[-1]).replace('t','g')
						if key in gene_id:
							print line
				
			

def call_gff(file,size):
	infile = file
	for chr in sorted(size.keys()):
	
		### store gene_ids
		gene_id = {}
		
		### hash coverage file
		hash={}
		print_gene_models = False
		key = COVERAGE
		hash_cov,gene_id = parse(file,chr,hash,key,print_gene_models,gene_id)
		
		### hash Augustus
		hash={}
		print_gene_models = False
		key = AUGUSTUS
		hash_aug,gene_id = parse(file,chr,hash,key,print_gene_models,gene_id)	
		

		### hash 454
		hash={}
		print_gene_models = False
		key = Key_454
		hash_TC,gene_id = parse(file,chr,hash,key,print_gene_models,gene_id)
		
		
		### hash the cufflinks gff co-ordinates
		hash={}
		key = CUFFLINKS
		print_gene_models = True
		hash,gene_id = parse(file,chr,hash,key,print_gene_models,gene_id) 
		filter_out(chr,hash,file)		
		
		
		### de-novo
		### hash the Augustus co-ordinates and find overlap against coverage file and TC
		print_gene_models = True
		file = 'temp'
		key = AUGUSTUS
		hash = {}
		hash,gene_id = find_overlap_denovo(file,chr,hash_cov,hash_TC,key,print_gene_models,gene_id)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		
		### hash the GeneMark co-ordinates and find overlap against coverage file and TC
		print_gene_models = True
		file = 'temp'
		key = GENEMARK
		hash = {}
		hash,gene_id = find_overlap_denovo(file,chr,hash_cov,hash_TC,key,print_gene_models,gene_id)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		
		### hash the GlimmerHMM co-ordinates and find overlap against coverage file and TC
		print_gene_models = True
		file = 'temp'
		key = GlimmerHMM
		hash = {}
		hash,gene_id = find_overlap_denovo(file,chr,hash_cov,hash_TC,key,print_gene_models,gene_id)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		'''
		### Homology based
		#os.system('touch protein_match.gff3')
		print_gene_models = True
		file = 'temp'
		key = AUGUSTUS
		hash = {}
		hash,gene_id = homology_evidence(file,chr,hash,print_gene_models,gene_id,key)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		
		### Homology based
		#os.system('touch protein_match.gff3')
		print_gene_models = True
		file = 'temp'
		key = GENEMARK
		hash = {}
		hash,gene_id = homology_evidence(file,chr,hash,print_gene_models,gene_id,key)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
		### Homology based
		#os.system('touch protein_match.gff3')
		print_gene_models = True
		file = 'temp'
		key = GlimmerHMM
		hash = {}
		hash,gene_id = homology_evidence(file,chr,hash,print_gene_models,gene_id,key)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		'''
		
		### hash the 454 co-ordinates and find overlap against coverage file and Augustus
		print_gene_models = True
		file = 'temp'
		key = Key_454
		hash = {}
		hash,gene_id = find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models,gene_id)
		os.system('cp temp temp2')
		file = 'temp2'
		filter_out(chr,hash,file)
		
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
	
	database = '../refseq.aa.fa'
	ref_seq = '../454Scaffolds.fna'
	
	#os.system('makeblastdb -in '+database +' >blast.temp')
	#os.system('formatdb -i '+ref_seq+' -p F -o T'+' >blast.temp')
	
	### get the end point for each chromosome
	size = get_size(file)
	
	### call gff file one by one
	call_gff(file, size)
	
	
	