'''
Created on Feb 25, 2013

@author: vgupta
'''

### script was made for filtering FASTQ file based on the read length range
### check output with 
# python ~/script/python/101_filter_fastq_len.py -i temp.fq -s 3 -l 5

import os,sys,getopt, re, commands
from C_loadFasta import *


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
    infile = ''
    min_len = 19
    max_len = 24
    
    try:
        opts, args = getopt.getopt(argv,"hi:s:l:",["ifile=","min_len","max_len"])
    except getopt.GetoptError:
        print '''
        python 101_filter_fastq_len.py \
                -i <inputfile> [multiple files separated by commas, no spaces]\
                -s <min_len> [19]\
                -l <max_len> [24]\
                -h <fastq_header> [@HWI] ## unique identifier for fastq reads
                
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
            python 101_filter_fastq_len.py \
                -i <inputfile> [multiple files separated by commas, no spaces]\
                -s <min_len> [19]\
                -l <max_len> [24]\
                -h <fastq_header> [@HWI] ## unique identifier for fastq reads 
            '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-s", "--min_len"):
            min_len = int(arg)
        elif opt in ("-l", "--max_len"):
            max_len = int(arg)
            
    logfile(infile)
    
    return infile, min_len, max_len

def filterFastq(file,min_len,max_len):
    count = 0
    read_sizes = {}
    ### open output file
    out=open(file+'.'+str(min_len)+'_'+str(max_len),'w')
    for line in open(file,'r'):
        if len(line) > 1:
            if line[0] != '#':
                line = line.strip()
                
                if line.startswith('@HWI'):
                    count = 0
                
                ### hash the reads
                count += 1
                if count%4 == 1:
                    read = line
                if count%4 == 2:
                    read += '\n'+line
                    size = len(line)
                    if size in read_sizes:
                        read_sizes[size] += 1
                    else:
                        read_sizes[size] = 1
                if count%4 == 3:
                    read += '\n'+line
                if count%4 == 0:
                    read += '\n'+line
                    
                ### print the read if it is in length range
                if count%4 == 0:
                    if min_len <= size <= max_len:
                        out.write(read+'\n')
    
    out.close()
    return read_sizes

if __name__ == "__main__":
    
    files,min_len,max_len = options(sys.argv[1:])
    
    for file in list(files.split(',')):
        ### filter fastq file
        file = file.strip()
        read_sizes = filterFastq(file,min_len,max_len)
        
        ### write number of lines in the fastq file
        o.write('Number of reads in the fastq file'+'\t'+file+'\t')
        o.write(str(int((commands.getoutput('wc -l '+file)).split()[0])/4)+'\n')
        
        ### write the length distribution
        o.write('Size Fractionation')
        for key in sorted(read_sizes):
            o.write(str(key)+'\t'+str(read_sizes[key])+'\n')
            
    
    ### close the logfile
    o.close()