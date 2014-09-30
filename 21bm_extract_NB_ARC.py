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
            python 100b_fasta2flat.py -i <infile>
            '''
    sys.exit(2)

def temp(file):
    return

### main argument to 

def options(argv):
    global infile, threads, lst
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:l:",["infile=","threads=","lst="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-l", "--lst"):
            lst = arg
            
    
    logfile(infile)
            

def LOADfasta(file):
    global headers, seq
    first_line = True
    seq = {}
    headers = {}
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
                headers[header] = ''
            else:
                string += line
            first_line = False
    if string != '':
        seq[header] = string
    return seq

def print_domain(lst):
    '''
    Q#1 - >AT1G61180.1[AT1G61180.1]	specific	250236	156	440	5.81933e-129	391.702	pfam00931	NB-ARC	 - 	cl18944
    '''
    last_header = ''
    for line in open(lst,'r'):
        line = line.strip()
        tokens = line.split('\t')
        if re.search('NB-ARC',line):
            header = line.split('[')[1].split(']')[0]
            start = int(tokens[3])
            end = int(tokens[4])
            if header in headers:
                if header != last_header:
                    print '>'+header
                    print seq[header][start-1:end]
                last_header = header

if __name__ == "__main__":
    

    options(sys.argv[1:])
        
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    LOADfasta(infile)
    
    print_domain(lst)
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)
