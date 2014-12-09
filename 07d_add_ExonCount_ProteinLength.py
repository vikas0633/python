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
    global infile, threads, occi, palle, repens
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:o:p:r:",["infile=","threads=","occi=","palle=","repens="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-o", "--occi"):
            occi = arg
        elif opt in ("-p", "--palle"):
            palle = arg
        elif opt in ("-r", "--repens"):
            repens = arg
    
    logfile(infile)
            
def hash_file(file):
    hash_temp = {}
    for line in open(file, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        if tokens[18] not in hash_temp:
            hash_temp[tokens[18]] = tokens[12]+"\t"+tokens[13]+"\t"+tokens[15]
        else:
            if int(tokens[12]) > int(hash_temp[tokens[18]].split("\t")[0]):
                hash_temp[tokens[18]] = tokens[12]+"\t"+tokens[13]+"\t"+tokens[15]
    return hash_temp

def process_file(occi_hash, palle_hash, repens_hash):
    for line in open(infile, 'r'):
        line = line.strip()
        tokens = line.split('\t')
        data1 = ''
        data2 = ''
        if tokens[1] in repens_hash:
            data1 = repens_hash[tokens[1]]
        if tokens[2] in occi_hash or tokens[2] in palle_hash:
            if tokens[2].startswith("occi"):    
                data2 = occi_hash[tokens[2].replace("occidentale","clover")]
            if tokens[2].startswith("palle"):
                data2 = palle_hash[tokens[2]]
        print tokens[0] +"\t"+tokens[1]+ "\t" + data1 + "\t" + data2

if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    occi_hash = hash_file(occi)
    palle_hash = hash_file(palle)
    repens_hash = hash_file(repens)
    
    process_file(occi_hash, palle_hash, repens_hash)
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

        

    