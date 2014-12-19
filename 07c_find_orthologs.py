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
    global infile, threads, blast, order
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:b:o:",["infile=","threads=","blast=","order="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
        elif opt in ("-b", "--blast"):
            blast = arg
        elif opt in ("-o", "--order"):
            order = arg
    
    logfile(infile)
            
def hash_blast():
    HashBlast = {}
    HashBlastPair = {}
    for line in open(blast,'r'): 
        line = line.strip()
        tokens = line.split("\t")
        key = tokens[0]+"_"+tokens[1]
        if key not in HashBlast:
            HashBlast[key] = "\t".join(tokens[2:])
        
        if tokens[0] not in HashBlastPair:
            HashBlastPair[tokens[0]] = [tokens[1]]
        else:
            HashBlastPair[tokens[0]].append(tokens[1])
            
    return HashBlast, HashBlastPair

def hash_contig_order():
    HashContig = {}
    for line in open(order, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        token = tokens[9].split('_')
        for i in range(len(token)/2):
            j = token[2*i] + '_' + token[2*i+1] 
            if tokens[0] not in HashContig:
                HashContig[tokens[0]] = ["g".join(j.split('g')[:-1])]
            else:
                HashContig[tokens[0]].append("g".join(j.split('g')[:-1]))
    
    return HashContig
        

def find_orthologs(HashBlast, HashBlastPair, HashContig):
    o=open(infile+'.nonOrtho','w')
    store_key1 = {}
    store_key2 = {}
    count = 0
    
    ### Condition -1: present in right contig and also has a blast hit
    for line in open(infile, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        contig_1 = "g".join(tokens[1].split('g')[:-1])
        contig_2 = "g".join(tokens[2].split('g')[:-1])    
        if len(tokens) >= 4:
            if tokens[3] != "No_Hit":  
                ### Not best reciprocal but still present in correct contig
                if tokens[1] not in store_key1 and tokens[2] not in store_key2:
                    if contig_1 in HashContig:
                        if int(tokens[3]) < 100:
                            for i in HashContig[contig_1]:
                                if contig_2 == i:
                                    if tokens[1]+"_"+tokens[2] in HashBlast:
                                        print "Condition-1:"+ '\t' + tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+HashBlast[tokens[1]+"_"+tokens[2]]
                                    else:
                                        print "Condition-1:"+ '\t' +tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+"Low quality blast"
                                    store_key1[tokens[1]] = ''
                                    store_key2[tokens[2]] = ''
                                    break
                            
    ### Condition-2: Best reciprocal hit
    for line in open(infile, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        contig_1 = "g".join(tokens[1].split('g')[:-1])
        contig_2 = "g".join(tokens[2].split('g')[:-1])
        if tokens[1] not in store_key1 and tokens[2] not in store_key2:
            if len(tokens) >= 4:
                if tokens[3] != "No_Hit":    
                    ### Best Reciprocal hit
                    if int(tokens[3]) == 1:
                        count += 1
                        print >> sys.stderr, count
                        if tokens[1]+"_"+tokens[2] in HashBlast:
                            print "Condition-2:"+ '\t' + tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+HashBlast[tokens[1]+"_"+tokens[2]]
                        else:
                            print "Condition-2:"+ '\t' +tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+"Low quality blast"
                        store_key1[tokens[1]] = ''
                        store_key2[tokens[2]] = ''
                        
    ### condition-3: find the best blast hit not previously used 
    for line in open(infile, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        contig_1 = "g".join(tokens[1].split('g')[:-1])
        contig_2 = "g".join(tokens[2].split('g')[:-1])    
        if len(tokens) >= 4:
            if tokens[3] != "No_Hit": 
                ### find the best blast hit not previously used
                if tokens[1] not in store_key1 and tokens[2] not in store_key2:
                    if tokens[1]+"_"+tokens[2] in HashBlast:
                        identity = float(HashBlast[tokens[1]+"_"+tokens[2]].split("\t")[0])
                        align_length = int(HashBlast[tokens[1]+"_"+tokens[2]].split("\t")[1])
                        #print identity, align_length
                        if identity > 90 or align_length>100:   
                            print "Condition-3:"+ '\t' +tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+HashBlast[tokens[1]+"_"+tokens[2]]
                            store_key1[tokens[1]] = ''
                            store_key2[tokens[2]] = ''
                            
    for line in open(infile, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        contig_1 = "g".join(tokens[1].split('g')[:-1])
        contig_2 = "g".join(tokens[2].split('g')[:-1])    
        ### find the best blast hit not previously used
        if tokens[1] not in store_key1:
            ### check of next blast hits
            if tokens[1] in HashBlastPair:
                for key in HashBlastPair[tokens[1]]:
                    if key not in store_key2:
                        if tokens[1]+"_"+key in HashBlast:
                            identity = float(HashBlast[tokens[1]+"_"+key].split("\t")[0])
                            align_length = int(HashBlast[tokens[1]+"_"+key].split("\t")[1])
                            #print identity, align_length
                            if identity > 90 or align_length>100:
                                print "Condition-4:"+ '\t' +tokens[1]+"\t"+key+"\t"+"0"+"\t"+HashBlast[tokens[1]+"_"+key]
                                store_key1[tokens[1]] = ''
                                store_key2[key] = ''
                                break
                                
    for line in open(infile, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        contig_1 = "g".join(tokens[1].split('g')[:-1])
        contig_2 = "g".join(tokens[2].split('g')[:-1])    
        if len(tokens) >= 4:
            if tokens[3] != "No_Hit": 
                ### find the best blast hit not previously used
                if tokens[1] not in store_key1 and tokens[2] not in store_key2:
                    if tokens[1]+"_"+tokens[2] in HashBlast:
                        identity = float(HashBlast[tokens[1]+"_"+tokens[2]].split("\t")[0])
                        align_length = int(HashBlast[tokens[1]+"_"+tokens[2]].split("\t")[1])
                        o.write("Condition-Fail:"+ '\t' +tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+HashBlast[tokens[1]+"_"+tokens[2]]+'\n')
                        store_key1[tokens[1]] = ''
                        store_key2[tokens[2]] = ''
                    else:
                        o.write("Condition-Fail:"+ '\t' +tokens[1]+"\t"+tokens[2]+"\t"+tokens[3]+"\t"+"No reverse Hit"+'\n')
                        store_key1[tokens[1]] = ''
                        store_key2[tokens[2]] = ''
        else:
            o.write("Condition-Fail:"+ '\t' +tokens[1]+"\t"+tokens[2]+"\t"+"No forward Hit"+'\n')
            store_key1[tokens[1]] = ''
            store_key2[tokens[2]] = ''
    o.close()
    o=open(infile+'.noBlastHit','w')
    for line in open(infile, 'r'):
        line = line.strip()
        tokens = line.split("\t")
        contig_1 = "g".join(tokens[1].split('g')[:-1])
        if tokens[1] not in store_key1:
            o.write("Condition-NoBlastHit:"+ '\t' +tokens[1]+"\t"+"No forward Hit"+'\n')
            store_key1[tokens[1]] = ''
    o.close()
if __name__ == "__main__":
    

    options(sys.argv[1:])
    
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    
    HashBlast, HashBlastPair = hash_blast()
    HashContig = hash_contig_order()
    
    find_orthologs(HashBlast, HashBlastPair, HashContig)
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)

    