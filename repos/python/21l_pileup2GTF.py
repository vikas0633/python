### 21l_pileup2GTF.py - script converts a pileup to a gtf file based on the coverage - /Users/vikas0633/Desktop/script/python




### options
'''
min_coverage cutoff - 5
min_fragment_length - 30
name as average coverage of the fragment
call island if fragment is "min_fragment_length":100 base pairs and there is nothing +/- 1 KB region
'''

### pileup file
'''
Ljchr1_pseudomol_20120830	425	C	1	^S.	b
Ljchr1_pseudomol_20120830	426	G	1	.	b
Ljchr1_pseudomol_20120830	427	A	1	.	b
'''

### gtf gile
'''
381 Twinscan  CDS          380   401   .   +   0  gene_id "001"; transcript_id "001.1";
381 Twinscan  CDS          501   650   .   +   2  gene_id "001"; transcript_id "001.1";
381 Twinscan  CDS          700   707   .   +   2  gene_id "001"; transcript_id "001.1";
381 Twinscan  start_codon  380   382   .   +   0  gene_id "001"; transcript_id "001.1";
381 Twinscan  stop_codon   708   710   .   +   0  gene_id "001"; transcript_id "001.1";
'''


import os, sys, getopt

def options(argv):
    infile = ''; min_cov = 5; min_len = 30; min_intron = 10;
    try:
    	opts, args = getopt.getopt(argv,"hp:c:l:i:",["pileup=","coverage=","min_length=","min_intron_length="])
    except getopt.GetoptError:
    	print 'python 21l_pileup2GTF.py -p <pileup> -c <coverage> -l <min_length_frag> -i <min_intron_length>'
    	sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print 'python 21l_pileup2GTF.py -p <pileup> -c <coverage> -l <min_length_frag> -i <min_intron_length>'
    		sys.exit()
    	elif opt in ("-p", "--pileup"):
    		infile = arg
    	elif opt in ("-c", "--coverage"):
    		min_cov = arg
    	elif opt in ("-l", "--min_length_frag"):
    		min_len = arg
    	elif opt in ("-i", "--min_intron_length"):
    		min_intron = arg
    return infile,min_cov, min_len, min_intron

def process_line(token,chr,start,end,coverage,pos,last_end,frag_no):
	
	### check if base qualifies the coverage criteria
	if int(token[3]) > int(min_cov): 
	
		### check start of new block
		if int(token[1]) - end > int(min_intron):
			if (end - start) > int(min_len):
				coverage = round(float(coverage)/pos,2)
				name = "exon"
				if  min_len < (end - start) < 100:
					if int(token[1]) - end > flank:
						if (end - last_end) > flank:
							name = "island"
				print chr+'\t'+'cufflinks'+'\t'+'exon'+'\t'+str(start)+'\t'+str(end)+'\t.\t+\t0\t'+'gene_id "cov_'+str(coverage)+'_'+name+'_'+str(frag_no)+'"; transcript_id cov_"'+str(coverage)+'_'+name+'_'+str(frag_no)+'.1";'
				last_end = end
				frag_no += 1
			coverage = 0
			start = int(token[1])
			pos = 0 
		pos += 1
		coverage += int(token[3])
		end = int(token[1])
		
	return start,end,coverage,pos,last_end, frag_no	
	


def make_gff():
	
	start = 0
	last_chr = ''
	first_line = True
	pos = 0
	last_end = start
	frag_no = 0
	### read the infile
	for line in open(infile,'r'):
		if line[0] != '#':
			line = line.strip()
			token = line.split('\t')
			
			if first_line == True: ### check for first chromosome
				chr = token[0]
				start = int(token[1])
				end = int(token[1])
				coverage = 0
				pos = 0
				last_end = start
				
			if token[0] != last_chr: ### check for new chromosome
				chr = token[0]
				start = int(token[1])
				end = int(token[1])
				coverage = 0
				pos = 0
				last_end = start
				
			start,end,coverage,pos,last_end,frag_no = process_line(token,chr,start,end,coverage,pos,last_end,frag_no)
			
			first_line = False		
			last_chr = chr


if __name__ == "__main__":
    
    infile,min_cov, min_len, min_intron = options(sys.argv[1:])  
    
    flank = 1000 ### flanking region for defining island
    make_gff()