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
    global infile, threads, gaps
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:g:",["infile=","threads=","gaps="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-g", "--gaps"):
            gaps = int(arg)
    
    logfile(infile)
            
def merge_contigs():
    first_line = True
    block_start = True
    for line in open(infile,'r'):
            
        line = line.strip()
        tokens = line.split('\t')
        contig = tokens[0]
        st = int(tokens[1])
        en = int(tokens[2])
        g_id = tokens[3].replace('ID=','')
        score = tokens[4]
        strand = tokens[5]
        
        if first_line == False:
            if contig == last_contig:
                ### check if new contig overlaps
                if bl_en + gaps > st: 
                    if en > bl_en:
                        bl_en = en
                        bl_id += '_' + g_id
                else:
                    block_start = True
                    print last_contig + '\t' + str(bl_st) + '\t' + str(bl_en) + '\t' + 'id=' + bl_id+ '\t' + str(bl_score) + '\t' + str(bl_strand)
            else:
                block_start = True
                print last_contig + '\t' + str(bl_st) + '\t' + str(bl_en) + '\t' + 'id=' + bl_id + '\t' + str(bl_score) + '\t' + str(bl_strand)
            
        ### check if new block starts
        if block_start == True:
            bl_st = st
            bl_en = en
            bl_id = g_id
            bl_contig = contig
            bl_score = score
            bl_strand = strand 
            block_start = False
            
        first_line = False
        last_contig = contig
        last_en = en
    
    ### print the last line
    print last_contig + '\t' + str(bl_st) + '\t' + str(bl_en) + '\t' + 'id=' + bl_id+ '\t' + str(bl_score) + '\t' + str(bl_strand)
        

if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Input file must be sorted"
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    merge_contigs()
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    