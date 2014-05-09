#-----------------------------------------------------------+
#                                                           |
# 21bg_find_fragmented_genemodels.py - Script to find fragmented genemodels               |
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


print 'Input GFF3 file must be sorted '

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
    global infile, max_dist, max_size, threads
    infile = ''
    max_dist = 1000
    max_size = 600
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
        elif opt in ("-m", "--max_dist"):
            max_dist = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
            
    
    logfile(infile)
            

def find_avg_gene_density(chromosome,count,gene_count_hash):
    correct_chro = False
    gene_count = 0
    for line in open(infile,'r'):
        line = line.strip()
        if len(line)>0 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            
            count += 1
            if count%10000 == 0:
                print 'Number of lines processed: ', chromosome, '{:9,.0f}'.format(count)
            
            if obj.seqids() == chromosome:
                correct_chro = True
                if obj.types()=='gene':
                    gene_count += 1
            if correct_chro == True and obj.seqids() != chromosome:
                break
            
    gene_count_hash[chromosome] = gene_count
    
def find_high_gene_density(chromosome,count,avg_gd):
    o_frag = open(infile+'.'+chromosome+'.frags.temp','w')
    last_start = -100000
    new_block = True
    correct_chro = False
    for line in open(infile,'r'):
        line = line.strip()
        if len(line)>0 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            
            if obj.seqids() == chromosome:
                correct_chro = True
                if obj.types()=='gene':
                    start = int(obj.starts())
                    
                    if start - last_start < max_dist:
                        if max_size > int(obj.ends()) - int(obj.starts()):
                            region += int(obj.ends()) - int(obj.starts())
                            gene_count += 1
                            gene_id.append(str(obj))
                            new_block = False
                    else:
                        
                        if new_block == False:
                            if gene_count*1000/float(gene_count) >= 5*avg_gd and gene_count>2:
                                print chromosome, block_start
                                for g_id in gene_id:
                                    o_frag.write(g_id+'\n')
                        
                        gene_count = 0
                        region = int(obj.ends()) - int(obj.starts())
                        new_block = True
                        gene_id = []
                        if max_size > int(obj.ends()) - int(obj.starts()):
                            gene_id.append(str(obj))
                            gene_count += 1
                        block_start = int(obj.starts())
                    
                    last_start = start
            if correct_chro == True and obj.seqids() != chromosome:
                break
    o_frag.close()

if __name__ == "__main__":
    

    file = options(sys.argv[1:])
    
    print 'Hashing the chromosomes name'
    chroHash = get_size(infile)
    
    print 'Chromosome sizes based on the GFF3 file are: '
    for chro in sorted(chroHash):
        print chro, chroHash[chro]
        
        
    ### multithreading        
    thread_list = []
    count = 0
    if len(chroHash) <= threads:
        manager = Manager()
        gene_count_hash = manager.dict()
        for chromosome in sorted(chroHash):
            count += 1
            t = Process(target=find_avg_gene_density, args=(chromosome,count,gene_count_hash))
            thread_list.append(t)
            t.start()
    
        for thread in thread_list:
            thread.join()
    else:
        gene_count_hash = {}
        for chromosome in sorted(chroHash):
            find_avg_gene_density(chromosome,count,gene_count_hash)
    
    print 'gene_count_hash:', gene_count_hash
    
    genome_size = sum(chroHash.values())
    gene_count_hash = dict(gene_count_hash)
    for i in sorted(gene_count_hash):
        print "Gene Density per kb: ", i, round(gene_count_hash[i]*1000/float(chroHash[i]),4)
    
    avg_gd = round(sum(gene_count_hash.values())*1000/float(genome_size),4)
    
    print 'Overall Gene density: ', avg_gd
    
    print 'Looking for the regions with more than 5 times gene density in 10 Kb windows'
    if len(chroHash) <= threads:
        for chromosome in sorted(chroHash):
            count += 1
            t = Process(target=find_high_gene_density, args=(chromosome,count,avg_gd))
            thread_list.append(t)
            t.start()
    
        for thread in thread_list:
            thread.join()
    else:
        for chromosome in sorted(chroHash):
            find_high_gene_density(chromosome,count,avg_gd)
    
    os.system('cat *.frags.temp > '+ infile+'.frags')
    os.system('rm *.frags.temp')
    
    ### close the logfile
    o.close()
    