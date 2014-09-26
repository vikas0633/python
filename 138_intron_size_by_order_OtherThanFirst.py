

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
            python 100b_fasta2flat.py -i <ifile>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global infile, threads
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:m:t:",["ifile=","max_dist=","threads="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
            
    
    logfile(infile)

## get PARENT ID
def get_PARENT(line):
    line = line.strip()
    match = re.search(r'Parent=.+',line)
    if match:
        return match.group().split(';')[0].replace('Parent=','')
    else:
        print 'Error at line'
        print line
        sys.exit('Parent ID is missing in the attributes')        

    
def print_intron(file):
    
    first_gene = True
    
    last_parent_ID = ''
    last_end = ''
    
    transcripts_pos = []
    
    for line in open(file,'r'):
        line = line.strip()
        if len(line) > 1 and not line.startswith('#'):
            tokens = line.split('\t')
            if tokens[2]!= "chromosome" and tokens[2]!= "protein"  and tokens[2]!= "five_prime_UTR" :
                obj_type = tokens[2]
                obj_strands = tokens[6]
                obj_start = int(tokens[3])
                obj_end = int(tokens[4])
                if obj_type == 'CDS': 
                    transcripts_pos.append(obj_start)
                    transcripts_pos.append(obj_end)
                if obj_type == 'mRNA' :
                    if len(transcripts_pos) > 2:
                        if first_gene == False:
                            transcripts_pos.sort()
                            if last_exon_strand == '+':
                                for i in range(1,len(transcripts_pos)/2 -1):
                                    print str(transcripts_pos[(2*i+2)] - transcripts_pos[(2*i+1)]) + '\t' + str(i+1)
                            if last_exon_strand == '-':
                                for i in range(len(transcripts_pos)/2 -2):
                                    print str(transcripts_pos[(2*i+2)] - transcripts_pos[(2*i+1)]) + '\t' + str(-1*i-1)
                    transcripts_pos = []
                    first_gene = False 
                    last_exon_strand = obj_strands     
        
    ###for last gene
    if len(transcripts_pos) > 2:
        if last_exon_strand == '+':
            for i in range(1,len(transcripts_pos)/2 -1):
                print str(transcripts_pos[(2*i+2)] - transcripts_pos[(2*i+1)]) + '\t' + str(i+1)
        if last_exon_strand == '-':
            for i in range(len(transcripts_pos)/2 -2):
                print str(transcripts_pos[(2*i+2)] - transcripts_pos[(2*i+1)]) + '\t' + str(-1*i-1)

            
    

if __name__ == "__main__":
    

    options(sys.argv[1:])

    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print_intron(infile)
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    