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

import sys,os, commands, getopt, re
import datetime
now = datetime.datetime.now()

global infile, prior

#o = open(str(now.strftime("%Y-%m-%d_%H%M_"))+'gene.gff','w')

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 100b_fasta2flat.py -i <ifile>
            '''
    sys.exit(2)

def options(argv):
	global infile, prior,database, genome, cores
	infile = ''
	prior = ''
	cores = 4
	try:
		opts, args = getopt.getopt(argv,"hi:p:d:g:c:",["ifile=","priority=","database=","genome=","cores="])
	except getopt.GetoptError:
		help()
	for opt, arg in opts:
		if opt == '-h':
			help()
		elif opt in ("-i", "--ifile"):
			infile = arg
		elif opt in ("-p", "--priority"):
			prior = arg
		elif opt in ("-d", "--database"):
			database = arg
		elif opt in ("-g", "--genome"):
			genome = arg
		elif opt in ("-c", "--cores"):
			cores = arg
	logfile(infile)
	return infile


def translate_dna(sequence):
 
    gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'', 'TAG':'',
    'TGC':'C', 'TGT':'C', 'TGA':'', 'TGG':'W',
    }
    proteinseq = ''
    for n in range(0,len(sequence),3):
        if gencode.has_key(sequence[n:n+3]) == True:
			proteinseq += gencode[sequence[n:n+3]]
    return proteinseq

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
						if (token[1] == 'CUFFLINKS') or  (token[1] == 'GlimmerHMM') or (token[1] == 'Augustus') or (token[1] == 'GENEMARK'):
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
		if len(line) > 0 and not line.startswith('#') and len(token)>3:
			if (token[1] == key) & (token[0] == chr) & (token[2] == "gene"):
				for i in range(int(token[3]),int(token[4])+1):
					temp[i]=''
				temp_set = set(temp)
				#if (len(temp_set.intersection(hash_cov_set)) > 0.8*(len(temp_set)))  or (len(temp_set.intersection(hash_aug_set)) > 0.8*(len(temp_set))): ### removing coverage condition
				if True:
					if (token[1] == 'GlimmerHMM') or (token[1] == 'Augustus') or (token[1] == 'GeneMark.hmm'):
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID]=''
						print line
					
					else:
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID]=''
						mRNA_id = ID.replace('path','mrna')
						gene_id[mRNA_id]=''
						print line
						
					
					for i in range(int(token[3]),int(token[4])+1):
						hash[i]=''
				
				elif token[1] == 'CUFFLINKS': ### let cufflinks gene models pass even if low coverage or no augustus
					ID=(line.split('=')[1]).split(';')[0]
					gene_id[ID]=''
					print line
					for i in range(int(token[3]),int(token[4])+1):
						hash[i]=''
					
	return hash,gene_id
		
def find_overlap_denovo(file,chr,hash_cov,key,print_gene_models,gene_id):
	hash = {}
	hash_cov_set = set(hash_cov)
	for line in open(file,'r'):
		temp = {}
		line = line.strip()
		token = line.split('\t')
		if len(line) > 1 and not line.startswith('#'):
			if (token[1] == key) & (token[0] == chr) & (token[2] == "gene"):
				for i in range(int(token[3]),int(token[4])+1):
					temp[i]=''
				temp_set = set(temp)
				if (len(temp_set.intersection(hash_cov_set)) > 0.8*(len(temp_set))): 
					if (token[1] == 'CUFFLINKS'):
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID]=''
						print line
					elif token[1] == "GlimmerHMM":
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID] = ''
						ID_tokens = ID.split('.')
						ID = 'model'+'.'+'.'.join(ID_tokens[1:])
						gene_id[ID]=''
						print line
					elif token[1] == "Augustus":
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID] = ''
						ID_tokens = ID.split('.')
						ID = 'model'+'.'+'.'.join(ID_tokens[1:])
						gene_id[ID]=''
						l = line.replace('gene.g','gene.g_'+chr+'_')
						print l
					elif (token[1] == "GeneMark.hmm"):
						ID=(line.split('=')[1]).split(';')[0]
						gene_id[ID] = ''
						ID = ID.replace('_g','_t')
						gene_id[ID]=''
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
	os.system('nice -n 19 blastp -num_threads '+cores+' -db '+database+' -query temp.fa -evalue 0.00000001 -outfmt 6 -out temp.blastout')
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
		
		if token[1] == "Augustus":
			ID=(line.split('=')[1]).split(';')[0]
			gene_id[ID] = ''
			ID_tokens = ID.split('.')
			ID = 'model'+'.'+'.'.join(ID_tokens[1:])
			gene_id[ID]=''
		
		l = line.replace('gene.g','gene.g_'+chr+'_')
		print l+'.homology'
		print >> sys.stderr, "Blasted:"+ID
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
							dna_temp = open('dna.temp.fa','w')
							temp_seq = ''
							for line2 in open('temp.fa','r'):
								line2 = line2.strip()
								if len(line2) > 1:
									if line2.startswith('>'):
										dna_temp.write(line2+'\n')
									else:
										temp_seq += line2 
							dna_temp.write(translate_dna(temp_seq)+'\n')
							dna_temp.close()
							os.system('cp dna.temp.fa temp.fa')
							print >> sys.stderr, "Blasting: "+line 
							hash,gene_id=do_blast(file,chr,hash,print_gene_models,gene_id,query)
						os.system('echo ">"'+str('temp')+' >temp.fa')
					new_gene = True
					query = line
				if token[2]=="exon":
					os.system('nice -n 19 fastacmd -d '+genome+' -p F -s '+token[0]+' -L '+str(token[3])+','+str(token[4]) +' | awk "NR>1" >>temp.fa')
			
	
	return hash,gene_id
	
	

def gene_stru(gene_id,infile,chr):
	
	### get Parent
	def get_parent(line):
	    match = re.search(r'Parent=.+',line)
	    if match:
	        return match.group().split(';')[0].replace('Parent=','')
	
	### process infile
	for line in open(infile,'r'):
		line = line.strip()
		token = line.split('\t')
		if len(line) > 0 and not line.startswith('#'):
			if (token[0]==chr):
				if (token[2] == 'mRNA') or (token[2] == 'exon'):
					if (token[1] == 'CUFFLINKS'):
						keys = (line.split('=')[-1]).split('.')
						key = keys[0]+'.'+keys[1]
						if key in gene_id:
							print line
					elif (token[1] == 'Augustus'):
						key = get_parent(line)
						if key in gene_id:
							l = line.replace('gene.g','gene.g_'+chr+'_')
							l = l.replace('model.g','model.g_'+chr+'_')
							print l
					elif (token[1] == 'GlimmerHMM'):
						match = re.search(r'Parent=.+',line)
						key = match.group().replace('Parent=','')
						if key in gene_id:
							print line
					elif (token[1] == 'GeneMark.hmm'):
						match = re.search(r'Parent=.+',line)
						key = match.group().replace('Parent=','')
						if key in gene_id:
							print line
					else:
						key = get_parent(line)
						if key in gene_id:
							print line
					
			

def call_gff(file,size):
	
	
	### get the priority order
	list_file_names = []
	hash_priors = {}
	first_line = True
	for line in open(prior, 'r'):
		line = line.strip()
		if first_line == True:
			names = line.split(',')
			for name in names:
				list_file_names.append(name)
			first_line = False
		else:
			line = line.strip()
			priors = line.split(',')
			for p in priors:
				hash_priors[p] = list_file_names[priors.index(p)]
			break
	
	infile = file
	for chr in sorted(size.keys()):
		
		print >> sys.stderr, chr
		
		### store gene_ids
		gene_id = {}
		
		### hash coverage file
		hash={}
		print_gene_models = False
		key = 'cufflinks'
		hash_cov,gene_id = parse(file,chr,hash,key,print_gene_models,gene_id)
		
		### hash Augustus ??Do I need to hash Augustus 
		hash={}
		print_gene_models = False
		key = 'Augustus'
		hash_aug,gene_id = parse(file,chr,hash,key,print_gene_models,gene_id)	
		
		
		for i in range(len(hash_priors)):
			evidence = hash_priors[str(i+1)].strip()
			
			
			
			
			if i != 0:
				file = 'temp'
			
			
			if evidence == 'RNA-seq':
				### check the evidence type
				### hash the cufflinks gff co-ordinates
				hash={}
				key = 'CUFFLINKS'
				print_gene_models = True
				hash,gene_id = find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models,gene_id) 
				filter_out(chr,hash,file)
			
			elif evidence == 'Augustus':
				### de-novo
				### hash the Augustus co-ordinates and find overlap against coverage file and TC
				print_gene_models = True
				key = 'Augustus'
				hash = {}
				hash,gene_id = find_overlap_denovo(file,chr,hash_cov,key,print_gene_models,gene_id)
				os.system('cp '+file+' temp2')
				file = 'temp2'
				filter_out(chr,hash,file)
	
				
			elif evidence == 'GeneMark':
				### hash the GeneMark co-ordinates and find overlap against coverage file and TC
				print_gene_models = True
				key = 'GeneMark.hmm'
				hash = {}
				hash,gene_id = find_overlap_denovo(file,chr,hash_cov,key,print_gene_models,gene_id)
				os.system('cp '+file+' temp2')
				file = 'temp2'
				filter_out(chr,hash,file)
		
			elif evidence == "Glimmer":
				### hash the GlimmerHMM co-ordinates and find overlap against coverage file and TC
				print_gene_models = True
				key = "GlimmerHMM"
				hash = {}
				hash,gene_id = find_overlap_denovo(file,chr,hash_cov,key,print_gene_models,gene_id)
				os.system('cp '+file+' temp2')
				file = 'temp2'
				filter_out(chr,hash,file)
			
			### if is a GMAP contigs
			else:
				print_gene_models = True
				key = str(i+1)
				hash = {}
				hash,gene_id = find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models,gene_id)
				os.system('cp '+file+' temp2')
				file = 'temp2'
				filter_out(chr,hash,file)
				
			'''
			if evidence == 'Augustus':
				### de-novo
				### hash the Augustus co-ordinates and find overlap against coverage file and TC
				print_gene_models = True
				key = 'Augustus'
				hash = {}
				hash,gene_id = homology_evidence(file,chr,hash,print_gene_models,gene_id,key)
				os.system('cp '+file+' temp2')
				file = 'temp2'
				filter_out(chr,hash,file)
			'''
			
			'''
				### Homology based
				#os.system('touch protein_match.gff3')
				print_gene_models = True
				file = 'temp'
				key = AUGUSTUS
				hash = {}
				hash,gene_id = homology_evidence(file,chr,hash,print_gene_models,gene_id,key)
				os.system('cp '+file+' temp2')
				file = 'temp2'
				filter_out(chr,hash,file)
				
				
				### Homology based
				#os.system('touch protein_match.gff3')
				print_gene_models = True
				file = 'temp'
				key = GENEMARK
				hash = {}
				hash,gene_id = homology_evidence(file,chr,hash,print_gene_models,gene_id,key)
				os.system('cp '+file+' temp2')
				file = 'temp2'
				filter_out(chr,hash,file)
				
				### Homology based
				#os.system('touch protein_match.gff3')
				print_gene_models = True
				file = 'temp'
				key = GlimmerHMM
				hash = {}
				hash,gene_id = homology_evidence(file,chr,hash,print_gene_models,gene_id,key)
				os.system('cp '+file+' temp2')
				file = 'temp2'
				filter_out(chr,hash,file)
				
				### hash the 454 co-ordinates and find overlap against coverage file and Augustus
				print_gene_models = True
				file = 'temp'
				key = Key_454
				hash = {}
				hash,gene_id = find_overlap(file,chr,hash_cov,hash_aug,key,print_gene_models,gene_id)
				os.system('cp '+file+' temp2')
				file = 'temp2'
				filter_out(chr,hash,file)
			'''
		
		### print gene structure
		gene_stru(gene_id,infile,chr)
		
		file = infile
		
		
def get_size(file):
	size = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 1:
			if line[0] != '#':
				token = line.split('\t')
				
				if token[0] not in size:
					size[token[0]] = int(token[4])
					
				else:
					if int(token[4]) > size[token[0]]:
						size[token[0]] = int(token[4])
				
	return size

if __name__ == "__main__":
	
	options(sys.argv[1:])
		
	#os.system('makeblastdb -in '+database +' >blast.temp')
	#os.system('formatdb -i '+ref_seq+' -p F -o T'+' >blast.temp')
	
	### get the end point for each chromosome
	size = get_size(infile)
	
	### call gff file one by one
	call_gff(infile, size)
	
	
	