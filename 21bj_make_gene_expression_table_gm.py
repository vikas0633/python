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
    return

### main argument to 

def options(argv):
    global infile, threads, corr
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:c:",["infile=","threads=","corr="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-c", "--corr"):
            corr = arg
            
    
    logfile(infile)
            

def Load_Expression():
    global header
    data = {}
    first_line = True
    for line in open(infile, 'r'):
        if first_line == True:
            header = line.strip()
        else:    
            line = line.strip()
            tokens = line.split('\t')
            data[tokens[0]] = '\t'.join(tokens[1:])
        first_line = False
    return data

def print_expression(data):
    global header
    gene_ID = {}
    first_line = True
    for line in open(corr,'r'):
        if first_line == False:
            line = line.strip()
            tokens = line.split('\t')
            try:
                if tokens[0] in data:
                    gene_ID[tokens[1]] = tokens[0]
            except:
                continue
        first_line = False 
    
    print "Gene_Id\tProbe_Id\t"+header
    for key in sorted(gene_ID):
        try:
            if not re.search('BLASTN',key) and not re.search('Ambiguous',key):
                print key+'\t'+gene_ID[key]+'\t'+data[gene_ID[key]]
        except:
            continue

if __name__ == "__main__":

    options(sys.argv[1:])
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)
    
    ### load the expression data
    data = Load_Expression()
    
    print_expression(data)