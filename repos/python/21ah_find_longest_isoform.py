'''
Created on March 06, 2013

@author: vgupta
'''

### script was made for finding longest isoform in the spider protein set
### 21ah_find_longest_isoform.py - script was made for finding longest isoform in the spider protein set - /Users/vikas0633/Desktop/script/python


import os,sys,getopt, re
from C_loadFasta import *


### write logfile

def logfile(infile):
    import datetime
    now = datetime.datetime.now()
    o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')
    o.write("Program used: \t\t%s" % "21ah_find_longest_isoform.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
    o.close()
            
    

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print 'python 21ah_find_longest_isoform.py -i <inputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python 21ah_find_longest_isoform.py -i <inputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
            
    return infile
    
def print_isoforms(headers,seqs,header_order): ### run through keys and print the fasta formatted file
    for key in header_order:
        if key in seqs:
            print headers[key]
            print seqs[key] 
            del seqs[key]
            del headers[key]
     

def longest_isoform(file):  ### find the longest isoform
    header_order = []
    first_line = True
    headers = {}
    seqs = {}
    for line in open(file,'r'):
        line = line.strip() 
        if line[0] == '>':          ### check for headers
            if first_line == False:
                if header in headers:
                    if len(seqs[header]) < len(seq):
                        headers[header] = header_name
                        seqs[header] = seq
                else:
                    headers[header] = header_name
                    seqs[header] = seq
                    
            seq = ''
            header = line[1:].split('_')[0]+'_'+line[1:].split('_')[3]### >L1_T1/1_Tarantula_WB_fr1,L1_T1/1_Tarantula_WB
            header_name = line
            header_order.append(header)
            first_line = False
        else:
            seq += line
            
    ### for last sequence
    if header in headers:
        if len(seqs[header]) < len(seq):
            headers[header] = header_name
            seqs[header] = seq
    else:
        headers[header] = header_name
        seqs[header] = seq
           
    print_isoforms(headers,seqs,header_order)

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    longest_isoform(file)