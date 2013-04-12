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

def sortGFF3(file,chr):
    hash = {}
    for line in open(file,'r'):
        line = line.strip()
        if len(line)>1:
            if line[0]!='#':
                token = line.split('\t')
                if token[0] == chr:
                    if token[2] == 'gene':                  ### get the gene name and use it start position as key
                        match = re.search(r'Name=.+',line)
                        if match:
                            match = match.group().split(';')[0].replace('Name=','')
                        else:
                            match = re.search(r'ID=.+',line)
                            match = match.group().split(';')[0].replace('ID=','')
                        match += token[3]
                        key = (int(token[3]),match)
                        hash[key] = line + '\n'
                    else:
                        hash[key] += line+'\n'
    for key in sorted(hash):
        print hash[key]

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    ### run it by chromosome
    size = E_get_chr_size_gff3.get_size(file)
    for chr in sorted(size):
        ### modify gene names
        sortGFF3(file,chr)
    ### close the logfile
    o.close()