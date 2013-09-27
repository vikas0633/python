#-----------------------------------------------------------+
#                                                           |
# 30b_AnnotationGenomicLocations.py                         |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
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
# python ~/script/python/30b_AnnotationGenomicLocations.py -i file_with_coords -g gff3_file


### import modules
import os,sys,getopt, re
import classGene
from C_loadFasta import *
import E_get_chr_size_gff3
import classGeneStructure

### global variables
global infile, chr_col, start_col, gff3, header

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "30b_AnnotationGenomicLocations.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 100b_fasta2flat.py   -i <ifile> [File with the positions to be annotated]
                                        -c 
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global infile, chr_col, start_col, gff3, header
    infile = ''
    gff3 = ''
    chr_col = 1
    start_col = 2
    header = True
    try:
        opts, args = getopt.getopt(argv,"hi:c:s:g:",["ifile=",'chr_col=','start_col=', 'gff3='])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-c", "--chr_col"):
            chr_col = int(arg)
        elif opt in ("-c", "--start_col"):
            start_col = int(arg)
        elif opt in ("-g", "--gff3"):
            gff3 = arg
            
    logfile(infile)

### split line
def split_line(line):
    return line.strip().split('\t')

### get ID
def get_ID(line):
    line = line.strip()
    match = re.search(r'ID=.+;',line)
    if match:
        return match.group().split(';')[0].replace('ID=','')
    else:
        print 'Error at line'
        print line
        sys.exit('ID is missing in the attributes')
    
### get PARENT ID
def get_PARENT(line):
    line = line.strip()
    match = re.search(r'Parent=.+',line)
    if match:
        return match.group().split(';')[0].replace('Parent=','')
    else:
        print 'Error at line'
        print line
        sys.exit('Parent ID is missing in the attributes')
        
def addCoords(gs, exons, cds, utr, intron, mRNA_type):
    
    exon_list = gs.GetExons()
    for i in exon_list:
        exons[i] = ''
    
    cds_list = gs.GetCDS()
    for i in cds_list:
        cds[i] = ''
    
    utr_list = gs.GetUTR()
    for i in utr_list:
        utr[i] = ''
        
    intron_list = gs.GetIntron()
    for i in intron_list:
        intron[i] = ''
    
    mRNA_type_list = gs.GetType()
    for i in mRNA_type_list:
        mRNA_type[i] = mRNA_type_list[i]
    
    return exons, cds, utr, intron, mRNA_type
                            
def hash_annotations(gff3, chro):
    exons = {}
    cds = {}
    utr = {}
    intron = {}
    inter = {}
    mRNA_type = {}
    first_line = True
    last_end = 0
    last_exon_end = 0
    parent_id = ''

    for line in open(gff3, 'r'):
        line = line.strip()
        if len(line) > 1 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            
            if obj.types() == 'mRNA':
                
                if first_line == False:
                    exons, cds, utr, intron, mRNA_type = addCoords(gs, exons, cds, utr, intron, mRNA_type)
                    
                first_line = False
                gs = classGeneStructure.GeneStructure(obj)
                gs.addmRNA(obj)
    
            if obj.types() == 'exon':
                gs.addexon(obj)
    
            if obj.types() == 'CDS':
                gs.addcds(obj)
    
    exons, cds, utr, intron, mRNA_type = addCoords(gs, exons, cds, utr, intron, mRNA_type)
    
    return exons, cds, utr, intron, mRNA_type

def annotate(infile, chro, exons, cds, utr, intron, mRNA_type):
    
    header = True
    for line in open(infile, 'r'):
        if len(line) > 1 and not line.startswith('#'):
            if header == False:
                line = line.strip()
                token = line.split('\t')
                chr = token[chr_col-1]
                start = int(token[start_col-1])
                if start in mRNA_type:
                    types = mRNA_type[start]
                else:
                    types = 'Intergenic'
                if chro == chr:
                    if start in cds:
                        print line+'\t'+ 'cds\t' + types
                    elif start in utr:
                        print line+'\t'+'utr\t' +  types
                    elif start in exons:
                        print line+'\t'+ 'exon\t' +  types
                    elif start in intron:
                        print line+'\t'+ 'intronic\t' +  types
                    else:
                        print line+'\t'+ 'intergenic\t' +  types
            else:
                header = False
        

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    chr_size = E_get_chr_size_gff3.get_size(gff3)
    
    for chro in sorted(chr_size):
        exons, cds, utr, intron, mRNA_type = hash_annotations(gff3, chro)
        
        ### annotate the points
        annotate(infile, chro, exons, cds, utr, intron, mRNA_type)
        
    ### close the logfile
    o.close()