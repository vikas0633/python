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
            python 141_trim_fasta.py -i <infile>
                -i/--infile STRING <input fasta>
                -t/--threads INT <number of cpus>
                -e/--end INT <5 or 3> "Operate on 5' or 3' end"
                -k/--keep INT <keep the bases> <default>
                -t/--trim INT <Remove these bases>
                
            '''
    sys.exit(2)

def temp(file):
    num_lines = sum(1 for line in open(file))
    return num_lines

### main argument to 

def options(argv):
    global infile, threads, end, keep, trim
    infile = ''
    threads = 2
    end = 5
    keep = 0
    trim = 0
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:e:k:t:",["infile=","threads=","end=","keep=","trim="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-e", "--end"):
            end = int(arg)
        elif opt in ("-k", "--keep"):
            keep = int(arg)
        elif opt in ("-t", "--trim"):
            trim = int(arg)
            
    
    logfile(infile)

def LOADfasta(file):
        o = open(file+'.temp','w')
	first_line = True
	seq = {}
	string = ''
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0 :			
			if line[0] == '>':
				if first_line == False:
					if string != '': 
						o.write('>'+header+'\n')
                                                o.write(string+'\n')
				string = ''
				header = line[1:].strip().split()[0]
			else:
				string += line
		first_line = False			
	if string != '':
                o.write('>'+header+'\n')
		o.write(string+'\n')
        o.close()
        
def process_fasta(seq):
    for line in open(infile+'.temp','r'):
        line = line.strip()
        if line.startswith('>'):
            print line
            header = line[1:].strip().split()[0]
        else:
            length = len(line)
            if end == 5:
                if keep != 0:
                    print line[:keep]
                elif trim != 0:
                    print line[trim:]
            elif end == 3:
                if keep != 0:
                    print line[length-keep:length]
                elif trim != 0:
                    print line[:length-trim]

if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    seq = LOADfasta(infile)
    
    process_fasta(seq)
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    