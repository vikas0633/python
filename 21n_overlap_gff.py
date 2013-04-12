### 21n_overlap_gff.py
### takes two or more gff files merge the files where you see an overlap


### python ~/Desktop/script/python/21n_overlap_gff.py sample.20121014_glimmerHMM_arabidopsis.gff3 sample.augustus.gff > overlap.gff

### sample gff
'''
Ljchr1_pseudomol_20120830	GlimmerHMM	mRNA	2404	3747	.	-	.	ID=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudomol_20120830.path1.gene1
Ljchr1_pseudomol_20120830	GlimmerHMM	CDS	2404	2514	.	-	0	ID=Ljchr1_pseudomol_20120830.cds1.1;Parent=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudomol_20120830.path1.gene1;Note=final-exon
Ljchr1_pseudomol_20120830	GlimmerHMM	CDS	2953	3108	.	-	0	ID=Ljchr1_pseudomol_20120830.cds1.2;Parent=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudomol_20120830.path1.gene1;Note=internal-exon
Ljchr1_pseudomol_20120830	GlimmerHMM	CDS	3198	3278	.	-	0	ID=Ljchr1_pseudomol_20120830.cds1.3;Parent=Ljchr1_pseudomol_20120830.path1.gene1;Name=Ljchr1_pseudom
'''

import sys

import datetime
now = datetime.datetime.now()

#o = open(str(now.strftime("%Y-%m-%d_%H%M_"))+'gene.gff','w')

def options(argv):
	
	gff = argv
	
	return gff



def parse(file,chr,hash):
	for line in open(file,'r'):
		line = line.strip()
		token = line.split('\t')
		if line[0] != '#':
			if (token[2] == 'mRNA') or (token[2] == 'gene'): 
				for i in range(int(token[3]),int(token[4])+1):
					hash.append(i)
	
	return hash

def make_gff(chr,pos):
	
	new_block = False
	start = -100
	end = -100
	
	for i in pos:
		
		if new_block == True:
			if last_start > 0:
				print chr+'\t'+'merged\t'+'gene\t'+str(last_start)+'\t'+str(last_end)+'\t.\t.\t.\t.'
			new_block = False
		if (i - end) > 10:
			new_block = True
			last_start = start
			last_end = end
			start = i
		end = i
		
def call_gff():
	for chr in size:
		### hash gff co-ordinates
		hash=[]
		for file in gff:
			hash = parse(file,chr,hash) 
		
		hash.sort()
		make_gff(chr,hash)
		
def get_size(file):
	size = {}
	for line in open(file,'r'):
		line = line.strip()
		if line[0] != '#':
			token = line.split('\t')
			
			if token[0] not in size:
				size[token[0]] = int(token[4])
				
			else:
				if int(token[4]) > size[token[0]]:
					size[token[0]] = int(token[4])
				
	return size

if __name__ == "__main__":
	
	gff = options(sys.argv[1:])  
	
	### get the end point for each chromosome
	size = get_size(gff[0])
	
	### call gff file one by one
	call_gff()
	
	