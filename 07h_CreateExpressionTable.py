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
    global infile, threads, ortho, occi, palle, repens
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:j:o:p:r:",["infile=","threads=","ortho=","occi=","palle=","repens="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-j", "--ortho"):
            ortho = arg
        elif opt in ("-o", "--occi"):
            occi = arg
        elif opt in ("-p", "--palle"):
            palle = arg
        elif opt in ("-r", "--repens"):
            repens = arg
            
    logfile(infile)
            
def hash_expression(FILE):
    hash_expression = {}
    first_line = True
    for line in open(FILE, 'r'):
        line = line.strip()
        if first_line == True:
            first_line = False
            header = line
        else:
            tokens = line.split()
            hash_expression[".".join(tokens[0].split('.')[:-1])] = '\t'.join(tokens[1:])
    return hash_expression, header

def hash_pairs(FILE):
    hash_pairs = {}
    first_line = True
    for line in open(FILE, 'r'):
        line = line.strip()
        tokens  = line.split('\t')
        hash_pairs[tokens[2]] = tokens[1]
    return hash_pairs

def print_out(infile, Occi_Hash, occi_header, Palle_Hash, palle_header, Repens_Hash, repens_header, Ortho):
    print occi_header + '\t' + repens_header + '\t' + palle_header + '\t' + '\t' + repens_header
    for line in open(infile, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        if tokens[0] in Ortho and tokens[1] in Ortho:
            print tokens[0] + '\t' + Occi_Hash[tokens[0]] + '\t' + Ortho[tokens[0]] + '\t' + Repens_Hash[Ortho[tokens[0]]] + '\t' + tokens[1] + '\t' + Palle_Hash[tokens[1]] + '\t' + Ortho[tokens[1]] + '\t' + Repens_Hash[Ortho[tokens[1]]]

if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    Occi_Hash, occi_header = hash_expression(occi)
    Palle_Hash, palle_header = hash_expression(palle)
    Repens_Hash, repens_header = hash_expression(repens)
    
    ### hash progenitor/hybrid pairs
    Ortho = hash_pairs(ortho)
    
    
    ### print out the expression of homeolougus pairs
    print_out(infile, Occi_Hash, occi_header, Palle_Hash, palle_header, Repens_Hash, repens_header, Ortho)
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    