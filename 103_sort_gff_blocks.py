'''
Created on Feb 25, 2013

@author: vgupta
'''


#103_sort_gff_blocks.py - script to sort GFF blocks - /Users/vikas0633/Desktop/script/python

import os,sys,getopt, re
from C_loadFasta import *
import E_get_chr_size_gff3


### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')

global count
count = 0

### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "103_sort_gff_blocks.py"+'\n')
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
            python 103_sort_gff_blocks.py -i <inputfile>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 103_sort_gff_blocks.py -i <inputfile>
                '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
            
    return infile

### split line
def split_line(line):
    return line.strip().split('\t')

### get ID
def get_ID(line):
    line = line.strip()
    match = re.search(r'ID=.+;',line)
    if match:
        return match.group().split(';')[0].replace('ID=','')
    else:
        print 'Error at line'
        print line
        sys.exit('ID is missing in the attributes')
    
### get PARENT ID
def get_PARENT(line):
    line = line.strip()
    match = re.search(r'Parent=.+;',line)
    if match:
        return match.group().split(';')[0].replace('Parent=','')
    else:
        print 'Error at line'
        print line
        sys.exit('ID is missing in the attributes')
    
def sortGFF3(file,chr):
    hash = {}
    gene_names = {} ### store  gene names
    gene_starts = {} ### store gene start
    gene = {} ### store gene line
    mRNA_names = {}
    global count
    for line in open(file,'r'):
        line = line.strip()
        if len(line)>1:
            if line[0]!='#':
                token = line.split('\t')
                if token[0] == chr:
                    if token[2] == 'gene':                  ### get the gene name and use it start position as key
                        match = get_ID(line)
                        gene_starts[int(token[3]), match] = match
                        gene[match] = line
                        gene_names[match] = ''
                        
                    elif token[2] == 'mRNA':
                        match = get_ID(line)
                        if match not in mRNA_names:
                            p_ID = get_PARENT(line)
                            if p_ID in hash:
                                hash[p_ID] += line + '\n'
                                
                            else:
                                hash[p_ID] = line + '\n'
                        
                        mRNA_names[match] = ''
                            
                    else:
                        hash[p_ID] += line + '\n'
                        
    
    for key in sorted(gene_starts):
        g_ID =  gene_starts[key]
        print gene[g_ID]
        print hash[g_ID]
    
if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    ### run it by chromosome
    size = E_get_chr_size_gff3.get_size(file)
    for chr in sorted(size):
        ### modify gene names
        sortGFF3(file,chr)
    ### close the logfile
    o.close()