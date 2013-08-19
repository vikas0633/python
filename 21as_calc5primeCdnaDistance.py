#-----------------------------------------------------------+
#                                                           |
# template.py - template for python scripting               |
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
from C_loadFasta import *
import E_get_chr_size_gff3


### global variables
global HEADER
HEADER='Chro\tStart\tEnd\tInsertID\tInsertCoord\tStrand\t5primeDistance\tTranscriptLength'
print HEADER
global infile, points

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
            python 100b_fasta2flat.py -i <ifile> [GFF3 file]
                                        -p <point> [genomic co-ordinates from Lore1, first 3 coulmns as chr, start, end]
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global infile, points
    infile = ''
    points = ''
    try:
        opts, args = getopt.getopt(argv,"hi:p:",["ifile=", "points="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-p", "--points"):
            points = arg
            
    logfile(infile)
            

def hash_coords_mRNA(coords_dis, exons, strand, coords_mRNA_len):
    
    ### process + strand
    ### from smallest exon to largest
    dist = 0
    if strand == '+':
        for i in range(len(exons)/2):
            for coords in range(exons[i], exons[i+1]+1, 1):
                dist += 1
                coords_dis[coords] = dist
        for i in range(len(exons)/2):
            for coords in range(exons[i], exons[i+1]+1, 1):
                coords_mRNA_len[coords] = dist
    
    ### from largest exon to smallest
    elif strand == '-':
        for i in range(len(exons)/2, 0 , -1):
            for coords in range(exons[i], exons[i-1]+1, -1):
                dist += 1
                coords_dis[coords] = dist
        for i in range(len(exons)/2, 0 , -1):
            for coords in range(exons[i], exons[i-1]+1, -1):
                coords_mRNA_len[coords] = dist
                
    return coords_dis, coords_mRNA_len

    
def hash_coords(file, chr):
    coords_dis={}
    coords_mRNA_len = {}
    exons = []
    first_gene = True
    for line in open(file, 'r'):
        if len(line) > 1 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            if obj.seqids() == chr:
                if obj.types() == "mRNA":
                    strand =  obj.strands()
                    if first_gene == False:
                        coords_dis, coords_mRNA_len = hash_coords_mRNA(coords_dis, exons, strand, coords_mRNA_len)
                        exons = []
                    first_gene = False
                if obj.types() == "exon":
                    exons.append(int(obj.starts()))
                    exons.append(int(obj.ends()))
            
    ### for last mRNA
    coords_dis, coords_mRNA_len = hash_coords_mRNA(coords_dis, exons, strand, coords_mRNA_len)
    
    return coords_dis, coords_mRNA_len

def get_5prime_dist(points, coords_dis, chr, coords_mRNA_len):
    
    for line in open(points):
        line = line.strip()
        if len(line) >1 and not line.startswith('#'):
            token = line.split('\t')
            if token[0] == chr:
                if int(token[1]) in coords_dis:
                    print line, coords_dis[int(token[1])], coords_mRNA_len[int(token[1])]


if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    chr_size = E_get_chr_size_gff3.get_size(infile)
    
    for chr in sorted(chr_size):
        ### Process file chromosomes
        
        ### hash all the genomic co-ordinates from gene's 5' end
        coords_dis, coords_mRNA_len = hash_coords(infile, chr)
        
        ### get the 5' distances
        get_5prime_dist(points, coords_dis, chr, coords_mRNA_len)
    
    ### close the logfile
    o.close()