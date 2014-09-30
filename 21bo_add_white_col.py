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
    global infile, threads, conf
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:c:",["infile=","threads=","conf="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-c", "--conf"):
            conf = arg
            
    
    logfile(infile)
            

def hash_colered():
    global gene_ID
    gene_ID = {}
    for line in open(conf,'r'):
        line = line.strip()
        tokens = line.split('\t')
        print line
        gene_ID[tokens[1]] = ''

def add_white():
    for line in open(infile,'r'):
        line = line.strip()
        if line.startswith('>'):
            line = line[1:]
            if re.search("Medicago",line):
                if len(line.split('|')) > 1:
                    g_id =  line.split('|')[1]
                    if g_id not in gene_ID:
                        print 'contain\t' + g_id + '\t' + str(255)+' '+str(255)+' '+str(255)+'\t'+'\t'+ '2'
            if re.search("Soybean",line):
                if len(line.split('|')) > 1:
                    g_id =  line.split('|')[0]
                    if g_id not in gene_ID:
                        print 'contain\t' + g_id + '\t' + str(255)+' '+str(255)+' '+str(255)+'\t'+'\t'+ '2'
            elif re.search(".",line):
                g_id =  line
                if g_id not in gene_ID:
                    print 'prefix\t' + g_id + '\t' + str(255)+' '+str(255)+' '+str(255)+'\t'+'\t'+ '2'
            else:
                g_id =  line
                if g_id not in gene_ID:
                    print 'prefix\t' + g_id + '\t' + str(255)+' '+str(255)+' '+str(255)+'\t'+'\t'+ '2'

if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    hash_colered()
    add_white()
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

        
    