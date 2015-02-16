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
            

def processBS(BS,CycleBS,line):
    cigar = line.split(':')[-1]
    cond = 'AEIPN'
    num = '' 
    for letter in cigar:
        if letter not in cond:
            num += letter
        else:
            if letter in CycleBS[BS]:
                CycleBS[BS][letter] += int(num)
            else:
                CycleBS[BS][letter] = int(num)
            num = ''
    return CycleBS
        
        
def count_CycleBS():
    count = 0; r_line=0; b_line=0; l_line=0;a_line=0;p_line=0;
    line1=''; line2='';line3='';line4='';
    CycleBS = {'#B':{},'#L':{},'#A':{},'#P':{}}
    for line in open(infile,'r'):
        line = line.strip()
        count += 1
        
        
        ### process ID line
        if count % 4 == 1:
            r_line += 1
            if re.search('#B',line):
                CycleBS = processBS('#B',CycleBS,line)
                b_line += 1
            elif re.search('#L',line):
                CycleBS = processBS('#L',CycleBS,line)
                l_line += 1
            elif re.search('#A',line):
                CycleBS = processBS('#A',CycleBS,line)
                a_line += 1
            elif re.search('#P',line):
                CycleBS = processBS('#P',CycleBS,line)
                p_line += 1
            
        ### process seq line
        
        ### process (strand ?) line

        ### process QUAL line        
        
    
    print "Total reads:\t", '{:9,.0f}'.format(r_line)  
    print "Live beads:\t", '{:9,.0f}'.format(l_line), CycleBS['#L']
    print "Blank reads:\t", '{:9,.0f}'.format(b_line), CycleBS['#B']
    print "Artifact reads:\t", '{:9,.0f}'.format(a_line), CycleBS['#A']
    print "Polyclonal reads:\t", '{:9,.0f}'.format(p_line), CycleBS['#P']
    
if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    #print >> sys.stderr, "Input count: " + str(temp(infile))
    
    count_CycleBS()
    
    #print >> sys.stderr, "Output count: " + str(temp(infile))
    #print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)


    