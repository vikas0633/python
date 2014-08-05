#-----------------------------------------------------------+
#                                                           |
# 31i_FixBoundries.py - script to modify GFF3 feature boundries                            |
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
            python 100b_fasta2flat.py -i <ifile>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global infile
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
    
    logfile(infile)
            
def FixBoundries():
    first_gene = True
    gene = []
    for line in open(infile,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            if obj.types() == 'gene':
                if first_gene == False:
                    obj = classGene.GFF3(gene[0])
                    print obj.seqids() + '\t' + \
                        obj.sources() + '\t' + \
                        obj.types() + '\t' + \
                        str(gene_start) + '\t' + \
                        str(gene_end) + '\t' + \
                        obj.scores() + '\t' + \
                        obj.strands() + '\t' + \
                        str(obj.phases()) + '\t' + \
                        obj.attributes()
                    for l in gene[1:]:
                        print l
                gene = []
                obj = classGene.GFF3(line)
                gene.append(line)
                gene_start = int(obj.starts())
                gene_end = int(obj.ends())
                first_gene = False
                
                
            elif obj.types() == 'mRNA':
                gene.append(line)
                if gene_start > int(obj.starts()):
                    gene_start = int(obj.starts())
                if int(obj.ends()) > gene_end:
                    gene_end = int(obj.ends())
                    
                mRNA_start = int(obj.starts())
                mRNA_end = int(obj.ends())
                mRNA_coord = {}
                for i in range(int(obj.starts()), int(obj.ends())+1):
                    mRNA_coord[i] = ''
            else:
                if int(obj.starts()) not in mRNA_coord:
                    line = obj.seqids() + '\t' + \
                    obj.sources() + '\t' + \
                    obj.types() + '\t' + \
                    str(int(mRNA_start)) + '\t' + \
                    obj.ends() + '\t' + \
                    obj.scores() + '\t' + \
                    obj.strands() + '\t' + \
                    str(obj.phases()) + '\t' + \
                    obj.attributes()
                if int(obj.ends()) not in mRNA_coord:
                    line = obj.seqids() + '\t' + \
                    obj.sources() + '\t' + \
                    obj.types() + '\t' + \
                    obj.starts() + '\t' + \
                    str(mRNA_end) + '\t' + \
                    obj.scores() + '\t' + \
                    obj.strands() + '\t' + \
                    str(obj.phases()) + '\t' + \
                    obj.attributes()
                        
                gene.append(line)
    
    obj = classGene.GFF3(gene[0])
    print obj.seqids() + '\t' + \
        obj.sources() + '\t' + \
        obj.types() + '\t' + \
        str(gene_start) + '\t' + \
        str(gene_end) + '\t' + \
        obj.scores() + '\t' + \
        obj.strands() + '\t' + \
        str(obj.phases()) + '\t' + \
        obj.attributes()
    
    for l in gene[1:]:
            print l
    
    
    
    
if __name__ == "__main__":
    

    options(sys.argv[1:])
    
 
    FixBoundries()
     
    o.close()
    