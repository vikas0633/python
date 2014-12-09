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
    global infile, threads, bam, contig_len
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:b:l:",["infile=","threads=","bam=","len="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-b", "--bam"):
            bam = arg
        elif opt in ("-l", "--len"):
            contig_len = arg
            
    
    logfile(infile)

def hash_st():
    HASH_ST = {}
    first_line = True
    for line in open(infile,'r'):
        line = line.strip()
        tokens = line.split()
        contig = tokens[0]
        st = int(tokens[1])
        en = int(tokens[2])
        g_id = tokens[3].replace('id=','')
        score = tokens[4]
        strand = tokens[5]
        if first_line == False and contig == last_contig:
            if last_g_id[0] != g_id:
                HASH_ST[last_contig, last_g_id, last_st] = st
        
        last_contig = contig
        last_st = st
        last_en = en
        last_g_id = g_id
        first_line = False
    return HASH_ST

def hash_map():
    HASH_MAP = {}
    for line in open(bam,'r'):
        line = line.strip()
        tokens = line.split('\t')
        g_id = tokens[0].replace('ID=','')
        if (tokens[1], g_id) not in HASH_MAP:
            HASH_MAP[tokens[1], g_id] = int(tokens[3])
        else:
            if HASH_MAP[tokens[1], g_id] < int(tokens[3]):
                HASH_MAP[tokens[1], g_id] = int(tokens[3])
    return HASH_MAP

def hash_len():
    HASH_LEN = {}
    for line in open(contig_len, 'r'):
        line = line.strip()
        tokens = line.split('\t')
        HASH_LEN[tokens[0]] = tokens[1]
    
    return HASH_LEN

def count_prog(HASH_ST, HASH_MAP, HASH_LEN):
    first_line = True
    first_model = True
    occi_bp = 0
    palle_bp = 0
    No_switch = 0
    Switch = 0
    Over_prog = 0
    Non_over_prog = 0
    last_contig = 'temp'
    last_g_id = 'temp'
    last_en = 1000000000000
    first_frag = True
    switches = ''
    print 'Contig\tContig_len\tOcci_bases\tPalle_bases\tOcci_fraction\tNo_switch\tSwitch\tOverlapping_progenators\tNon-Overlapping_progenators\tswitches'
    for line in open(infile,'r'):
        line = line.strip()
        tokens = line.split()
        contig = tokens[0]
        st = int(tokens[1])
        en = int(tokens[2])
        g_id = tokens[3].replace('id=','')
        score = tokens[4]
        strand = tokens[5]
            
        if first_line == False and contig != last_contig:
            if occi_bp > 0 or palle_bp > 0:
                print last_contig + '\t' + str(HASH_LEN[last_contig]) + '\t' + str(occi_bp) + '\t' + str(palle_bp) + '\t' + str(round(float(occi_bp)/(occi_bp+palle_bp),3)) + '\t' + str(No_switch) + '\t' + str(Switch) + '\t' + str(Over_prog) + '\t' + str(Non_over_prog) + '\t' + str(switches)
            occi_bp = 0
            palle_bp = 0
            No_switch = 0
            Switch = 0
            Over_prog = 0
            Non_over_prog = 0
            first_frag = True
            first_model = True
            switches = ''
            last_g_id = 'temp'
            last_en = 1000000000000
        if en - st > 1000:
            if re.search('occi',line):
                occi_bp += en - st
            if re.search('palle',line):
                palle_bp += en - st
        
            if contig == last_contig and first_frag == False:
                ### Count the switch/Non-switch
                #print last_g_id, g_id
                if (last_g_id[0:4] != g_id[0:4]) and (st > last_en) and (last_g_id != 'temp'):
                    if (contig, g_id, st) in HASH_ST:
                        if HASH_ST[contig, g_id, st] > en:
                            ### check for mapping quality
                            max_map = 0
                            gene_id_tokens = g_id.split('_')
                            for i in range(len(gene_id_tokens)/2):
                                gene_id = gene_id_tokens[2*i] + '_' + gene_id_tokens[2*i+1]
                            if HASH_MAP[contig, gene_id] > max_map:
                                max_map = HASH_MAP[contig, gene_id] 
                        
                                if max_map > 100:
                                    #print 'accepted: ' + last_g_id +'_'+ g_id
                                    Switch += 1
                                    last_en = en
                                    switches += last_g_id + '->' + g_id
                                    last_g_id = g_id
                if last_g_id[0] == g_id[0]:
                    No_switch += 1
                    last_en = en
                    last_g_id = g_id
                #print last_g_id[0] , g_id[0], No_switch, Switch  
                ### Count overlapping fragments
                if st < last_en:
                    Over_prog += 1
                else:
                    Non_over_prog += 1
            
            if first_frag == True:
                if (contig, g_id, st) in HASH_ST:
                    if HASH_ST[contig, g_id, st] > en:
                        ### check for mapping quality
                        max_map = 0
                        gene_id_tokens = g_id.split('_')
                        for i in range(len(gene_id_tokens)/2):
                            gene_id = gene_id_tokens[2*i] + '_' + gene_id_tokens[2*i+1]
                            if HASH_MAP[contig, gene_id] > max_map:
                                max_map = HASH_MAP[contig, gene_id] 
                                if max_map > 100:
                                    last_g_id = g_id
                                    first_frag = False
            
            last_contig = contig
        
        first_line = False
    if occi_bp > 0 or palle_bp > 0:
        print last_contig + '\t' + str(HASH_LEN[last_contig]) + '\t' + str(occi_bp) + '\t' + str(palle_bp) + '\t' + str(round(float(occi_bp)/(occi_bp+palle_bp),3)) + '\t' + str(No_switch) + '\t' + str(Switch) + '\t' + str(Over_prog) + '\t' + str(Non_over_prog) + '\t' + str(switches)
    
if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    HASH_ST = hash_st()
    HASH_MAP = hash_map()
    HASH_LEN = hash_len()
    
    count_prog(HASH_ST, HASH_MAP, HASH_LEN)
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    