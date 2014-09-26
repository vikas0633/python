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
    
def file_empty(file):
    count = sum([1 for line in open(file)])
    if count == 0:
        sys.exit(file+' is empty')
    else:
        return count

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
    
    logfile(infile)
            
def find_reciprocal(infile):
    
    count_lines = 0 
    ### hash the blast hit against target using score column
    hash_target = {}
    hash_target_score = {}
    hash_target_query = {}
    for line in open(infile,'r'):
        line = line.strip()
        tokens = line.split('\t')
        score = float(tokens[11])
        count_lines += 1
        if count_lines%100000 == 0:
                print >> sys.stderr, str('Hashing blast output: ' + '{:9,.0f}'.format(count_lines))
        
        ## store the score values
        if tokens[1] not in hash_target:
            hash_target[tokens[1]] = line
            hash_target_score[tokens[1]] = score
            hash_target_query[tokens[1]] = tokens[0]

        elif score > hash_target_score[tokens[1]]:
            hash_target[tokens[1]] = line
            hash_target_score[tokens[1]] = score
            hash_target_query[tokens[1]] = tokens[0]
            
    count_lines = 0
    count_121_hits = 0
    query=''
    for line in open(infile,'r'):
        line = line.strip()
        tokens = line.split('\t')
        count_lines += 1
        if count_lines%100000 == 0:
                print >> sys.stderr, str('Searching for reciprcals: ' + '{:9,.0f}'.format(count_lines))
        
        if tokens[0] in hash_target_query.values():
            if tokens[0] == hash_target_query[tokens[1]] and query != tokens[0]:
                count_121_hits += 1
                print line + '\t' +  hash_target[tokens[1]]
        
        ### make sure you only check first hit
        query = tokens[0]
    
    return count_121_hits

if __name__ == "__main__":
    options(sys.argv[1:])
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Start finding 1-to-1 blast hits: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input Blast-hit pairs: " + str(file_empty(infile))
    print >> sys.stderr, "Total 1-to-1 reciprcal pairs: " + str(find_reciprocal(infile))
    print >> sys.stderr, "Completed finding 1-to-1 blast hits: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)
    
 
    o.close()
    