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


def get_size(file):
    size = {}
    for line in open(file,'r'):
        line = line.strip()
        if len(line) > 0:
            if line[0] != '#':
                token = line.split('\t')
                
                if token[0] not in size:
                    size[token[0]] = int(token[4])
                    
                else:
                    if int(token[4]) > size[token[0]]:
                        size[token[0]] = int(token[4])
    size_sorted={}
    for w in sorted(size, key=size.get, reverse=False):
        size_sorted[w]=size[w]
    return size_sorted


### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 100b_fasta2flat.py -i <infile>
            '''
    sys.exit(2)

def temp(file):
    num_lines = sum(1 for line in open(file))
    return num_lines

### main argument to 

def options(argv):
    global cuff, threads, aug
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hc:t:a:",["cuff=","threads=","aug="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-c", "--cuff"):
            cuff = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-a", "--aug"):
            aug = arg
            
            

def hash_cuff():
    Hash_Cuff = {}
    for line in open(cuff,'r'):
        line = line.strip()
        tokens = line.split('\t')
        if len(tokens) > 3:
            if tokens[2] == 'gene':
                if tokens[0] in Hash_Cuff:
                    Hash_Cuff[tokens[0]] += 1
                else:
                    Hash_Cuff[tokens[0]] = 1
    return Hash_Cuff

def hash_aug():
    Hash_Aug = {}
    for line in open(aug,'r'):
        line = line.strip()
        tokens = line.split('\t')
        if len(tokens) > 3:
            if tokens[2] == 'gene':
                if tokens[0] in Hash_Aug:
                    Hash_Aug[tokens[0]] += 1
                else:
                    Hash_Aug[tokens[0]] = 1
    return Hash_Aug

def print_bact(Hash_Cuff,Hash_Aug):
    print 'Contig' + '\t' + 'Augustus genes' + '\t' + 'Cufflinks genes'
    for key in Hash_Aug:
        if key in Hash_Cuff:
            if Hash_Aug[key] > 5*Hash_Cuff[key] and Hash_Aug[key] > 50:
                print key + '\t' + str(Hash_Aug[key]) + '\t' + str(Hash_Cuff[key]) 

if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    
    Hash_Cuff = hash_cuff()
    Hash_Aug = hash_aug()
    
    print_bact(Hash_Cuff,Hash_Aug)
    
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    