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
    global infile, threads, fasta
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:f:",["infile=","threads=","fasta="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-f", "--fasta"):
            fasta = arg
            
    
    logfile(infile)

### C_loadFasta.py - script to load fasta sequences 
def HashFasta(fasta):
    first_line = True
    seq = {}
    string = ''
    for line in open(fasta,'r'):
            line = line.strip()
            if len(line) > 0 :			
                    if line[0] == '>':
                            if first_line == False:
                                    if string != '' and len(string) == 2001: 
                                            seq[header] = string
                            string = ''
                            header = '_'.join(line[1:].strip().split()[0].split('_')[0:6])
                    else:
                            string += line
            first_line = False			
    if string != '' and len(string) == 2001: 
            seq[header] = string
    return seq

def print_homoeo_pair(seq):
    os.system('rm -rf 01_homoeologs')
    os.system('mkdir -p 01_homoeologs') 
    count = 0
    for line in open(infile, 'r'):
        line = line.strip()
        tokens = line.split('\t')
        if len(tokens) == 2:
            if tokens[0] in seq and tokens[1] in seq:
                count += 1
                outfile = '/array/users/vgupta/14_cloverAnalysis/08_regele/01_homoeologs/'+str(str(count).zfill(4))+'_homoeolog'+'.fa'  
                o = open(outfile,'w')
                o.write('>'+tokens[0]+'\n')
                o.write(seq[tokens[0]]+'\n')
                o.write('>'+tokens[1]+'\n')
                o.write(seq[tokens[1]]+'\n')
                o.close()
                os.system('nice -n 19 clustalw -infile='+outfile+' -ALIGN -TYPE=DNA -OUTFILE='+outfile+'.aln '+'-STATS='+outfile+'.stat')

if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    print >> sys.stderr, "Output count: " + str(temp(infile))
    
    seq = HashFasta(fasta)
    
    print_homoeo_pair(seq)
    
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    