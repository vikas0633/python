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
    global infile, threads, coverage
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:c:",["infile=","threads=","coverage="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-c", "--coverage"):
            coverage = int(arg)
            
    
    logfile(infile)
            

def LOADfasta(file):
    headers = []
    first_line = True
    seq = {}
    string = ''
    for line in open(file,'r'):
            line = line.strip()
            if len(line) > 0 :			
                    if line[0] == '>':
                            if first_line == False:
                                    if string != '': 
                                            seq[header] = string
                            string = ''
                            header = line[1:].strip().split()[0]
                            headers.append(header)
                    else:
                            string += line
            first_line = False			
    if string != '': 
            seq[header] = string
    return seq, headers

def count_coverage(seqs, headers):
    
    ### count coverage of each position
    pos_cov = {}
    out = open(infile+'.'+str(coverage)+'.fa',"w")
    for seq in seqs:
        for i in range(len(seqs[seq])):
            if i in pos_cov:
                if seqs[seq][i]=='-':
                    continue
                else:
                    pos_cov[i] += 1
            else:
                if seqs[seq][i]=='-':
                    continue
                else:
                    pos_cov[i] = 1
    
    ### print positions with higher coverage
    for seq in headers:
        out.write('>'+seq+'\n')
        count = 0
        for i in range(len(seqs[seq])):
            if pos_cov[i] > coverage:
                count += 1
                out.write(seqs[seq][i])
                if count%60 == 0:
                    out.write('\n')
        out.write('\n')
    for i in range(len(pos_cov)):
        print i, pos_cov[i]
            


if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    print 'Hashing the chromosomes name'
    
    seqs, headers = LOADfasta(infile)
    count_coverage(seqs, headers)
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    