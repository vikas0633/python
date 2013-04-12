#21k_make_input4_glimmerHMM.py - this scripts takes a gene structure file (gff3) and makes a exon file parsable by glimmerHMM - /Users/vikas0633/Desktop/script/python

import os,sys, getopt


### python 21k_make_input4_glimmerHMM.py -p 20121011_culinks.orf -g 05_transcripts.gff3 -r /u/vgupta/lotus_3.0/Ljchr0-6_pseudomol_20120830.chlo.mito.fa
### gtf gile
'''
##gff-version 3
Ljchr2_pseudomol_20120830	Cufflinks	gene	42766555	42767521	.	.	.	ID=CUFF.40101;
Ljchr2_pseudomol_20120830	Cufflinks	transcript	42766555	42767521	.	.	.	ID=CUFF.40101.1;Parent=CUFF.40101;
Ljchr2_pseudomol_20120830	Cufflinks	exon	42766555	42767521	1000	.	.	ID=exon:CUFF.40101.1:1;Parent=CUFF.40101.1;
Ljchr2_pseudomol_20120830	Cufflinks	gene	23151253	23151792	.	.	.	ID=CUFF.37175;
Ljchr2_pseudomol_20120830	Cufflinks	transcript	23151253	23151792	.	.	.	ID=CUFF.37175.1;Parent=CUFF.37175;
Ljchr2_pseudomol_20120830	Cufflinks	exon	23151253	23151792	1000	.	.	ID=exon:CUFF.37175.1:1;Parent=CUFF.37175.1;'''
### orf file

'''
#ID	orf	Frame	Start	Codon	Stop	Codon	Length
Ljchr2_pseudomol_20120830:15743098-15744028	1	+2	296	ATG	889	TAG	594
Ljchr2_pseudomol_20120830:20361288-20362966	1	-2	1440	ATG	94	TAG	1347
Ljchr2_pseudomol_20120830:20225413-20227205	1	+2	113	ATG	757	TGA	645
Ljchr2_pseudomol_20120830:40139625-40143937	1	-1	2410	ATG	1361	TAG	1050
Ljchr2_pseudomol_20120830:33122344-33123767	1	-2	1299	ATG	604	TAG	696
'''

def options(argv):
    infile = ''
    try:
    	opts, args = getopt.getopt(argv,"hp:g:r:",["orf_file=","gene_structure_file=","ref_seq="])
    except getopt.GetoptError:
    	print 'python 21k_make_input4_glimmerHMM.py -p <orf_file> -g <gene_structure_file> -r <ref_seq>'
    	sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print 'python 21k_make_input4_glimmerHMM.py -p <orf_file> -g <gene_structure_file> -r <ref_seq>'
    		sys.exit()
    	elif opt in ("-p", "--orf_file"):
    		orf_file = arg
    	elif opt in ("-g", "--gene_structure_file"):
    		gene_structure_file = arg
    	elif opt in ("-r", "--ref_seq"):
    		ref_seq = arg
    return orf_file,gene_structure_file, ref_seq

def hash_orf(infile):
	orf = {}
	for line in open(infile,'r'):
		line = line.strip()
		token = line.split('\t')
		orf[token[0]]=''
	return orf
	
	
def print_glimmer(orf,file,header):
	flag = False
	for line in open(file,'r'):
		line = line.strip()
		token = line.split('\t')
		if len(token) > 3:
			if (token[2] == "gene"):
				key = token[0]+':'+str(int(token[3])-1)+'-'+token[4].strip()
				if key in orf:
					if (token[0] == header):
						print
						flag = True
				else:
					flag = False
			if (token[2] == "exon") & (flag == True):
				print header+'\t'+token[3]+'\t'+token[4]
	
		

def call_by_chro(ref_seq,orf,gene_structure_file):
	for line in open(ref_seq,'r'):
		line=line.strip()
		if line[0] == '>':
			header = line[1:]
			### print glimmerHMM friendly file
			print_glimmer(orf,gene_structure_file,header)

if __name__ == "__main__":
    
    orf_file,gene_structure_file,ref_seq = options(sys.argv[1:])    
    
    ### hash the orf file to get the gene positions
    orf = hash_orf(orf_file)
    
    ### process ref seq fasta file
    call_by_chro(ref_seq,orf,gene_structure_file) 