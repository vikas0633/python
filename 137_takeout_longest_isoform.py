#-----------------------------------------------------------+
#                                                           |
# template.py - Python template                             |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                        |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# Short script to convert and copy the wheat BACs           |
# Run this in the parent dir that the HEX* dirs exist       |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/100b_fasta2flat.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame

### import modules
import os,sys,getopt, re

import threading

from multiprocessing import Process, Queue, Manager
from threading import Thread
import classGene
### global variables
global infile

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')

def file_empty(file):
    count = sum([1 for line in open(file)])
    if count == 0:
        sys.exit(file+' is empty')
    else:
        return count

def seq_count(file):
    count = 0
    for line in open(file,'r'):
        if line.startswith('>'):
            count += 1
    
    return count


### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 100b_fasta2flat.py -i <ifile>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global infile, threads
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
            
    
def find_longest_isoform(infile):
    global temp_header, temp_seq, gene_seq, gene_seq_len
    def store_gene():
        if temp_header in gene_seq:
            if len(temp_seq) > gene_seq_len[temp_header]:
                gene_seq[temp_header] = temp_seq
                gene_seq_len[temp_header] = len(temp_seq)
        else:
            gene_seq[temp_header] = temp_seq
            gene_seq_len[temp_header] = len(temp_seq)

    first_line = True
    gene_seq = {}
    gene_seq_len = {}
    for line in open(infile,'r'):
        line = line.strip()
        
        if line.startswith('>'):
            
            if first_line == False:
                store_gene()
                
            temp_header = line[1:].split('.')[0]
            temp_seq = ''
            
        else:
            temp_seq += line
        
        first_line = False
    ### for last seq
    store_gene()
    for header in gene_seq:
        print '>'+header
        print gene_seq[header]
    
    return len(gene_seq)

if __name__ == "__main__":
    

    options(sys.argv[1:])
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Finding longest isoform: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input transcript count: " + str(seq_count(infile))
    print >> sys.stderr, "Output gene count: " + str(find_longest_isoform(infile))
    print >> sys.stderr, "Completed finding longest isoform: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)
     
    o.close()
    