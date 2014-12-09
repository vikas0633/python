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
import operator
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
    global infile, threads
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:",["infile=","threads="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
            
    
    logfile(infile)
            
def count_genes():
    count_contig = 1
    count_same_contig = 0
    count_different_contig = 0
    count_replace = 0
    count_replaced = 0
    new_contig = True
    first_line = True
    count = 1
    gene = {}
    temp = {}
    last_contig_1 = ''
    last_contig_2 = ''
    for line in open(infile,'r'):
        line  = line.strip()
        token = line.split('\t')
        if len(line) > 1:
            contig_1 = "g".join(token[0].split('g')[:-1])
            contig_2 = "g".join(token[1].split('g')[:-1])
            if len(token)>3:
                if last_contig_1 == contig_1:
                    count  += 1
                    if token[1] not in temp:
                        if contig_2 in gene:
                            gene[contig_2] += 1
                        else:
                            gene[contig_2] = 1
                    temp[token[1]] = ''
                else:
                    if first_line == False:
                        print last_contig_1 +"\t"+ max(gene.iteritems(), key=operator.itemgetter(1))[0] +"\t"+ str(gene[max(gene.iteritems(), key=operator.itemgetter(1))[0]])
                        gene = {}
                        temp = {}
                        count  = 1
                        gene[contig_2] = 1
                        temp[token[1]] = ''
                
                    last_contig_1 = contig_1
                    last_contig_2 = contig_2
        first_line = False
    print last_contig_1 +"\t"+ max(gene.iteritems(), key=operator.itemgetter(1))[0] + '\t' + str(gene[max(gene.iteritems(), key=operator.itemgetter(1))[0]])

if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    
    count_genes()
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    