## 21m_gff2genestru.py
## script creates input for gb format conversion script
## this script is made to convert genome based gene models to gene based gene models
## I know it sounds confusing but in the output each gene is considered as a seperate 
## sequence and then gene structure file will have respective exon-intron structure 
## with respect to gene co-ordinates




### genome file
'''
>Ljchr1_pseudomol_20120830
CGAAACCCTGAAACTCTAAATACCCGAAACCCTAAAGCTCTGAAACCATGAAAACCCTGA
TTCTCGAAACCTCGAAACGCTCGAACCCCGAAACCCTAACGCTATGAAACCCTAAAACTA
TAAATACCTGAAACCCTAAAGCTTCGAAACCATGAAAACCCTGATTCTCGAAACCCAGAA
ATGCTCGAACACTGAAACTCTAAATACCCGAAACGGTCGAACCCCGAAACCATGAAAACC

'''

### cufflinks gff file
'''
Ljchr1_pseudomol_20120830	Cufflinks	gene	1899	5503	.	.	.	ID=CUFF.25163;
Ljchr1_pseudomol_20120830	Cufflinks	transcript	1899	5503	.	.	.	ID=CUFF.25163.1;Parent=CUFF.25163;
Ljchr1_pseudomol_20120830	Cufflinks	exon	1899	5503	1000	.	.	ID=exon:CUFF.25163.1:1;Parent=CUFF.25163.1;
Ljchr1_pseudomol_20120830	Cufflinks	gene	2089	15721	.	-	.	ID=CUFF.25164;
Ljchr1_pseudomol_20120830	Cufflinks	transcript	2089	15721	.	-	.	ID=CUFF.25164.1;Parent=CUFF.25164;
'''


import os, sys, getopt

def options(argv):
    gff = ''; ref = 5; 
    try:
    	opts, args = getopt.getopt(argv,"hg:f:",["gene_model=","fasta="])
    except getopt.GetoptError:
    	print 'python 21l_pileup2GTF.py -g <gene_model> -f <fasta> '
    	sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print 'python 21l_pileup2GTF.py -g <gene_model> -f <fasta> '
    		sys.exit()
    	elif opt in ("-g", "--gene_model"):
    		gff = arg
    	elif opt in ("-f", "--fasta"):
    		ref = arg
    return gff,ref


def make_gff():
	### gene must be first feature of each gene structure
	
	import datetime
	now = datetime.datetime.now()
	
	#o = open(str(now.strftime("%Y-%m-%d_%H%M_"))+'gene.gff','w')
	o = open(str(now.strftime("CUFF."))+'gene.gff','w')
	for line in open(gff,'r'):
		if line[0] != '#':
			line = line.strip()
			token = line.split()
			
			if token[2] == 'gene':
				gene_id = line.split('=')[1][:-1]
				start = int(token[3])
				new_gene = True
				
			if new_gene == True:
				o.write(gene_id+'\t'+token[1]+'\t'+token[2]+'\t'+str(int(token[3])-start)+'\t'+str(int(token[4])-start)+'\t'+token[5]+'\t'+token[6]+'\t'+token[7]+'\t'+token[8]+'\n')
		
	o.close()
				
				
def make_fasta():
	import datetime
	now = datetime.datetime.now()
	
	### make database for ref seq
	os.system('formatdb -i '+ref+' -p F -o T')
	
	#file = str(now.strftime("%Y-%m-%d_%H%M_"))+'gene.fa'
	file = str(now.strftime("CUFF."))+'gene.fa'
	
	o = open(file,'w')
	
	o.close()
	for line in open(gff,'r'):
		line = line.strip()
		token = line.split()
		if line[0] != '#':
			if token[2] == 'gene':
				gene_id = '>'+line.split('=')[1][:-1]
				## write the header
				os.system("echo '"+gene_id+"' >> "+ file )
				### use fasta command to fetch the sequences
				os.system('fastacmd -d '+ref+' -p F -s '+token[0]+' -L '+token[3]+','+token[4]+'| sed 1d >> '+ file) 
	

if __name__ == "__main__":
	
	gff,ref = options(sys.argv[1:])  
	
	### gff file with respect to gene co-ordinates
	make_gff()
	
	### make fasta file with gene sequences
	### using cufflinks gene names
	make_fasta()