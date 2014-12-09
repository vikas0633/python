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
import operator
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
            
def count_genes():
    count_contig = 1
    count_oo = 0
    count_op = 0
    count_po = 0
    count_pp = 0
    last_contig_1 = ''
    last_contig_2 = ''
    new_contig = True
    first_line = True
    gene = {}
    for line in open(infile,'r'):
        line  = line.strip()
        token = line.split('\t')
        contig_1 = "g".join(token[0].split('g')[:-1])
        contig_2 = "g".join(token[1].split('g')[:-1])
        prog = contig_2.split('_')[0]
        if len(token)==2:
            if last_contig_1 == contig_1:
                count_contig += 1
                
                if last_prog == prog:
                    if last_prog == 'occidentale':
                        count_oo += 1
                    elif last_prog == 'pallescens':
                        count_pp += 1
                else:
                    if last_prog == 'occidentale':
                        count_op += 1
                    elif last_prog == 'pallescens':
                        count_po += 1
            
            last_contig_1 = contig_1
            last_contig_2 = contig_2
            last_prog = prog
            first_line = False

    print "Total number of contigs proccessed: \t"+ str(count_contig)
    print "Progenator: \tTO\tTP"
    print "TO\t" + str(count_oo) + '\t' + str(count_op)
    print "TP\t" + str(count_po) + '\t' + str(count_pp)
    
    count_contig = 1
    count_o = 0
    count_p = 0
    count_ow = 0
    count_pw = 0
    last_contig_1 = ''
    last_contig_2 = ''
    new_contig = True
    first_line = True
    gene = {}
    gene_count = {}
    gene_count_no = 1
    for line in open(infile,'r'):
        line  = line.strip()
        token = line.split('\t')
        contig_1 = "g".join(token[0].split('g')[:-1])
        contig_2 = "g".join(token[1].split('g')[:-1])
        prog = contig_2.split('_')[0]
        if len(token)==2:
            if last_contig_1 == contig_1:
                gene_count_no += 1
                count_contig += 1
                if contig_2 in gene:
                    gene[contig_2] += 1
                else:
                    gene[contig_2] = 1

            if last_contig_1 != contig_1 and first_line == False:
                if len(gene) > 1:    
                    gene_max = max(gene.iteritems(), key=operator.itemgetter(1))[0] ## return the contig with highest gene count
                    
                    for key in gene:
                        if key == gene_max:
                            if key.startswith('occidentale'):
                                count_o += gene[key]
                            elif key.startswith('pallescens'):
                                count_p += gene[key]
                        else:
                            if key.startswith('occidentale'):
                                count_ow += gene[key]
                            elif key.startswith('pallescens'):
                                count_pw += gene[key]
                if gene_count_no not in gene_count:
                    gene_count[gene_count_no] = 1
                else:
                    gene_count[gene_count_no] += 1
                gene = {}
                gene_count_no = 1
                if contig_2 in gene:
                    gene[contig_2] += 1
                else:
                    gene[contig_2] = 1
            last_contig_1 = contig_1
            last_contig_2 = contig_2
            last_prog = prog
            first_line = False    
    print "Total number of contigs proccessed: \t"+ str(count_contig)
    print "Genes in single contig in order from occidentale:\t" + str(count_o) 
    print "Genes in single contig in order from pallenscens:\t" + str(count_p) 
    print "Genes in single contig from worng order from occidentale:\t" + str(count_ow)    
    print "Genes in single contig from worng order from pallescens:\t" + str(count_pw)  

    print "Gene counts in the contigs"
    print gene_count
    
    ### wrong genes at the end
    count_contig = 1
    count_o = 0
    count_p = 0
    count_ow = 0
    count_pw = 0
    last_contig_1 = ''
    last_contig_2 = ''
    new_contig = True
    first_line = True
    gene = {}
    gene_count = []
    gene_count_no = 1
    for line in open(infile,'r'):
        line  = line.strip()
        token = line.split('\t')
        contig_1 = "g".join(token[0].split('g')[:-1])
        contig_2 = "g".join(token[1].split('g')[:-1])
        prog = contig_2.split('_')[0]
        if len(token)==2:
            if last_contig_1 == contig_1:
                count_contig += 1
                if contig_2 in gene:
                    gene[contig_2] += 1
                else:
                    gene[contig_2] = 1
                gene_count.append(contig_2)
            if last_contig_1 != contig_1 and first_line == False:
                if len(gene) > 1:    
                    gene_max = max(gene.iteritems(), key=operator.itemgetter(1))[0] ## return the contig with highest gene count
                    
                    key = gene_count[-1]
                    if key != gene_max:
                        if key.startswith('occidentale'):
                            count_o += 1
                        elif key.startswith('pallescens'):
                            count_p += 1
                    if key == gene_max:
                        if key.startswith('occidentale'):
                            count_ow += 1
                        elif key.startswith('pallescens'):
                            count_pw += 1
                gene_count = []
                gene_count.append(gene_count_no)
                gene = {}
                if contig_2 in gene:
                    gene[contig_2] += 1
                else:
                    gene[contig_2] = 1
                gene_count.append(contig_2)
            last_contig_1 = contig_1
            last_contig_2 = contig_2
            last_prog = prog
            first_line = False    
    print "Wrong gene at the edge from occidentale:\t" + str(count_o) 
    print "Wrong gene at the edge from pallescens:\t" + str(count_p) 
    print "Right gene at the edge from occidentale:\t" + str(count_ow)
    print "Right gene at the edge from occidentale:\t" + str(count_pw)  

if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    
    count_genes()
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    