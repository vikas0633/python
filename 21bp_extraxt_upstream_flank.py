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

### make a script to extract N-kb fragments upstream of a gene given gene_ID, genome fasta file and GFF3 file.

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
    global infile, threads, gff3, flank, db
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:f:g:d:",["infile=","threads=","flank=","gff3=","database="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-f", "--flank"):
            flank = int(arg)
        elif opt in ("-g", "--gff3"):
            gff3 = arg
        elif opt in ("-d", "--database"):
            db = arg
        
    
    logfile(infile)
    
### store upstream closest gene
def store_upstream_gene():
    last_chro = ''
    last_g_id = ''
    last_start = 0
    last_end = 0
    last_last_end = -100000000
    last_last_start = -100000000
    last_strand = ''
    hash_upstream_gene = {}
    for line in open(gff3, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        if tokens[2] == "gene":
            chro = tokens[0]
            start = int(tokens[3])
            end = int(tokens[4])
            strand = tokens[6]
            g_id = line.split('"')[1]
            if last_chro == chro:
                if last_strand == "+":
                    if flank > last_start - last_last_end and last_last_start < last_start:
                        if last_start - last_last_end < 30:
                            if flank - (last_start - last_last_start) < 0:
                                hash_upstream_gene[last_g_id] = [last_chro, last_start - flank, last_start - flank]
                            else:
                                hash_upstream_gene[last_g_id] = [last_chro, last_start - flank, last_last_start]
                        else:
                            hash_upstream_gene[last_g_id] = [last_chro, last_last_end, last_start ]
                    else:
                        if start - flank > 0:
                            hash_upstream_gene[last_g_id] = [last_chro, start - flank, start]
                        else:
                            hash_upstream_gene[last_g_id] = [last_chro, 0, start]
                elif last_strand == "-":
                    if flank > start - last_end and end > last_end:
                        if start - last_end < 30:
                            if last_end + flank - end < 0:
                                hash_upstream_gene[last_g_id] =  [last_chro,  last_end + flank, last_end + flank]
                            else:
                                hash_upstream_gene[last_g_id] = [last_chro,  end, last_end + flank]
                        else:
                            hash_upstream_gene[last_g_id] = [last_chro, last_end, start]
                    else:
                        hash_upstream_gene[last_g_id] =  [last_chro, end, end + flank]
            elif last_chro != chro:
                last_last_end = -100000000
                last_last_start = -100000000
                last_start = 0
                last_end = 0
            last_last_end = last_end
            last_last_start = last_start
            last_end = end
            last_start = start
            last_chro  = chro
            last_g_id = g_id
            last_strand = strand
    return hash_upstream_gene

### store GeneID, strand and start pos        
def store_gene():
    hash_gene = {}
    for line in open(gff3, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        if tokens[2] == "gene":
            g_id = line.split('"')[1]
            strand = tokens[6]
            chro = tokens[0]
            if strand == "+":
                start = int(tokens[3])
            elif strand == "-":
                start = int(tokens[4])
            hash_gene[g_id] = [chro, strand, start]
    return hash_gene

def process_genes(hash_gene):
    for line in open(infile,'r'):
        line = line.strip()
        g_id = line.split()[0]
        [chro, start, end] = hash_gene[g_id]
            
        os.system('echo ">"'+str(g_id)+"_"+str(chro)+"_"+str(start)+"_"+str(end)+' >>'+infile+".fa")
        os.system('nice -n 19 fastacmd -d '+db+' -p F -s '+chro+' -L '+str(start)+','+str(end) +'| awk "NR>1" >>'+infile+".fa")

        

if __name__ == "__main__":
    

    options(sys.argv[1:])
    

    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    hash_upstream_gene = store_upstream_gene()
    #hash_gene =  hash_upstream_gene()
    process_genes(hash_upstream_gene)
    
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    
    o.close()
    