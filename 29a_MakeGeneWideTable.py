
'''
Created on June 06, 2013

@author: vgupta
'''

### script was made for making genewide table of snpEff
# 29a_MakeGeneWideTable.py - /Users/vikas0633/Desktop/script/python/ - script to put the snpEff data togehter

# Usage 
# 

import os,sys,getopt, re
import sys
import os.path
import glob

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
    directory = ''
    suff= ''
    gene_ID = ''
    try:
        opts, args = getopt.getopt(argv,"hd:s:g:",["dir=","suff=","gene_ID="])
    except getopt.GetoptError:
        print '''
            python 29a_MakeGeneWideTable.py \
                                            -d <dir> \
                                            -s <suff> \
                                            -g <gene_ID>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
            python 29a_MakeGeneWideTable.py \
                                            -d <dir> \
                                            -s <suff> \
                                            -g <gene_ID>
                '''
            sys.exit()
        elif opt in ("-d", "--dir"):
            directory = arg
        elif opt in ("-s", "--suff"):
            suff = arg
        elif opt in ("-g", "--gene_ID"):
            gene_ID = arg

            
    logfile(directory)
            
    return directory, suff, gene_ID
    
                            
### Loop over the files in a folder
def LoopFiles(direcotry, suff):
    for root, dirs, files in os.walk(direcotry):
        for file in files:
            if file.endswith(suff):
                
                process(file)
                
                
### process the gene_ID information
def GeneID(gene_ID):
    gene_id_hash = {}
    for line in open(gene_ID,'r'):
        if len(line) > 1:
            if line[0] != '#':
                token = line.strip().split('\t')
                gene_id_hash[token[4]] = "\t".join(token[0:4])
    print len(gene_id_hash)
    return gene_id_hash

if __name__ == "__main__":
    direcotry, suff, gene_ID = options(sys.argv[1:])
    
    ### process the gene_ID information
    gene_id_hash = GeneID(gene_ID) 
    
    ### loop over the files in the dir
    LoopFiles(direcotry, suff)
    
    ### close the logfile
    o.close()
    
    

