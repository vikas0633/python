### 21i_RMoutput2GTF.py - take a tab-formatted RMoutput file as parse it to make a gtf file - /Users/vikas0633/Desktop/script/python

''''
GTF stands for Gene transfer format. It borrows from GFF, but has additional structure that warrants a separate definition and format name.
Structure is as GFF, so the fields are: 
<seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

Here is a simple example with 3 translated exons. Order of rows is not important.

381 Twinscan  CDS          380   401   .   +   0  gene_id "001"; transcript_id "001.1";
381 Twinscan  CDS          501   650   .   +   2  gene_id "001"; transcript_id "001.1";
381 Twinscan  CDS          700   707   .   +   2  gene_id "001"; transcript_id "001.1";
381 Twinscan  start_codon  380   382   .   +   0  gene_id "001"; transcript_id "001.1";
381 Twinscan  stop_codon   708   710   .   +   0  gene_id "001"; transcript_id "001.1";
The whitespace in this example is provided only for readability. In GTF, fields must be separated by a single TAB and no white space.
'''


import os,sys, getopt


def options(argv):
    infile = ''; tab = True
    try:
    	opts, args = getopt.getopt(argv,"hi:n:",["ifile=","not_tab_format="])
    except getopt.GetoptError:
    	print 'python 21i_RMoutput2GTF.py -i <inputfile>'
    	sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print 'python 21i_RMoutput2GTF.py -i <inputfile>'
    		sys.exit()
    	elif opt in ("-i", "--ifile"):
    		infile = arg
    	elif opt in ("-n", "--not_tab_format"):
    		tab = False  
    return infile, tab
    
def convert2gtf(file): 
	for line in open(file,'r'):
		if line[0] in '0123456789':
			line = line.strip()
			token = line.split('\t')
			print str(token[4])+'\t'+str("RepeatMasker")+'\t'+str(token[10])+'\t'+str(token[5])+'\t'+str(token[6])+'\t'+str('.')+'\t'+str('+')+'\t'+str('0')+'\t'+'gene_id "'+str(token[9])+'_'+str(token[14])+'"; '+'transcript_id "'+str(token[9])+'_'+str(token[14])+'.1";'

if __name__ == "__main__":
    
    infile,tab = options(sys.argv[1:])    
    
    ### parse infile to make gtf file
    convert2gtf(infile)