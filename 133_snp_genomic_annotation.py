#-----------------------------------------------------------+
#                                                           |
# 133_snp_genomic_annotation.py - script to calculate the snps genomic distribution               |
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
import classGene

from multiprocessing import Process, Queue, Manager
from threading import Thread
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
    global infile, gff3
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:g:",["ifile=","gff3="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-g", "--gff3"):
            gff3 = arg
            
    logfile(infile)

def hash_GFF3(chromosome):
    
    count = 0
    coords = {}
    for line in open(gff3,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            if obj.seqids() == chromosome:
                count += 1
                if count%1000 == 0:
                    print 'Number of lines processed: ', chromosome, '{:9,.0f}'.format(count)
                for i in range(int(obj.starts()), int(obj.ends())+1):
                    if obj.types()=='mRNA':
                        coords[i] = 'mRNA'
                    if obj.types()=='exon':
                        coords[i] = 'exon'
                    if obj.types()=='CDS':
                        coords[i] = 'CDS'                        
    return coords

    
def annotate_snp(chromosome, snp_count):
    
    coords = hash_GFF3(chromosome)
    
    genic = 0
    exonic = 0
    cds = 0
    utr = 0
    intronic = 0
    intergenic = 0
    
    count = 0
    for line in open(infile,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            token = line.split()
            if token[0] == chromosome:
                count += 1
                if count%1000 == 0:
                    print 'Number of lines processed: ', chromosome, '{:9,.0f}'.format(count)
                if int(token[1]) in coords:
                    genic += 1
                    
                    if coords[int(token[1])] == 'CDS':
                        cds += 1
                        exonic += 1
                        
                    elif coords[int(token[1])] == 'exon':
                        exonic += 1
                        utr += 1
                    
                    else:
                        intronic += 1
                    
                else:
                    intergenic += 1
            
    snp_count[chromosome] =  str(genic) + '\t' + str(exonic) + '\t' + str(cds) + '\t' + str(utr) + '\t' + str(intronic) + '\t' + str(intergenic)     
    
if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    print 'Hashing the chromosomes name'
    chroHash = get_size(gff3)
    
    ### multithreading
    thread_list = []
    manager = Manager()
    snp_count = manager.dict()
    for chromosome in sorted(chroHash):
        t = Process(target=annotate_snp, args=(chromosome,snp_count))
        t.start()
        thread_list.append(t)
    for thread in thread_list:
        thread.join()
    
    snp_count = dict(snp_count)
    
    ### print results
    print 'chromosome'+'\t'+'genic' + '\t' + 'exonic' + '\t' + 'cds' + '\t' + 'utr' + '\t' + 'intronic' + '\t' + 'intergenic'
    for chro in snp_count:
        print chro+'\t'+snp_count[chro]
    
    ### close the logfile
    o.close()