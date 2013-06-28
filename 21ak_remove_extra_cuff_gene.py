'''
Created on April 05, 2013

@author: vgupta
'''

### 21ak_remove_extra_cuff_gene.py - script to filter out extra gene lines in the GFF3 file due to alternative transcripts from Cufflinks - /Users/vikas0633/Desktop/script/python


import os,sys,getopt, re


### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    

### main argument to 

def options(argv):
    ''' str -> (list)
    helper function for parsing options
    '''
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print '''
            python 100b_fasta2flat.py -i <inputfile>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 100b_fasta2flat.py -i <inputfile>
                '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
    
    return infile
    
                            
def remove_genes(file):
    ''' file -> (NULL)
    this function removes unnecessary gene names assigned by cuyfflinks 
    '''   
    gene_ID_hash = {}
    for line in open(file,'r'):
        line = line.strip()
        if len(line) > 2:
            if line[0] != '#':
                token = line.split('\t')
                if token[2] == "gene":
                    gene = line
                    gene_ID = re.search(r'ID=.+',line)
                    gene_ID = gene_ID.group().split(';')[0].replace('ID=','') 
                    print gene_ID, last_gene_ID
                    if gene_ID not in last_gene_ID:
                        print line
                        last_gene_ID = gene_ID
                    
                else:
                    print line
                    

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    
    ### parse the file and remove extra genes
    remove_genes(file)
    
    ### close the logfile
    o.close()