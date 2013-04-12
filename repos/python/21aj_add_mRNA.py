### 21aj_add_mRNA.py - script to add dummy mRNAs if absent - /Users/vikas0633/Desktop/script/python

'''
Created on March 27, 2013

@author: vgupta
'''

## Input
## chr5    rRNA    gene    7457324    7457488    .    -    .    ID=chr5_7822648_7822812#+.path5;Name=chr5_7822648_7822812#+;Type=rRNA

import os,sys,getopt, re
from C_loadFasta import *



### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "21aj_add_mRNA.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print '''
            python 21aj_add_mRNA.py -i <inputfile>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 21aj_add_mRNA.py -i <inputfile>
                '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
            
    return infile

def print_mRNA(last_gene):
    last_gene_token = last_gene.split('\t')
    gene_ID = re.search(r'ID=.+',last_gene)
    gene_ID = gene_ID.group().split(';')[0].replace('ID=','') 
    gene_type = re.search(r'Type=.+',last_gene)
    gene_type = gene_type.group().split(';')[0].replace('Type=','')
    print '\t'.join(last_gene_token[0:2])+'\tmRNA\t'+'\t'.join(last_gene_token[3:8])+'\t'+'ID='+gene_ID+'.1;'+'Parent='+gene_ID+';'+'Type='+gene_type+';'

                            
### this is a fnction to add dummy mRNAs
def dummy_mRNA(file):
    first_gene = True
    for line in open(file,'r'):
        line = line.strip()
        if len(line) > 1:
            if line[0] != '#':
                token = line.split('\t')
                
                ### check if we hit a gene
                if token[2]=="gene":
                    if first_gene == False:
                        if mRNA_flag == True:               ### check if mRNA already has been printed  
                            print_mRNA(last_gene)
                    print line
                    mRNA_flag = True
                    first_gene = False
                    gene = line
                    last_gene = gene
                elif token[2]=="mRNA":
                    print line
                    mRNA_flag = False
                else:
                    print line
    if mRNA_flag == True:               ### check if mRNA already has been printed  
        print_mRNA(last_gene)
        
if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    ### add dummy to genes
    dummy_mRNA(file)
    
    ### close the logfile
    o.close()