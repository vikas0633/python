'''
Created on Feb 25, 2013

@author: vgupta
'''

### script was made for formatting the longest ORF to flat file

### check output with 
# python ~/script/python/100b_fasta2flat.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame

import os,sys,getopt, re
from C_loadFasta import *

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print 'python 100b_fasta2flat.py -i <inputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python 100b_fasta2flat.py -i <inputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    return infile
    
def makeflat(file):
    o = open('flat_file.txt','w')
    for line in open(file,'r'):
        line = line.strip()
        if line[0]=='>':
            o.write('\n'+line.split('_fr')[0][1:] +'\t')
        else:
            o.write(line)            
    o.close()

if __name__ == "__main__":
    file = options(sys.argv[1:])
    
    makeflat(file) 