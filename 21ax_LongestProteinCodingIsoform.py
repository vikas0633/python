#-----------------------------------------------------------+
#                                                           |
# 21ax_LongestProteinCodingIsoform.py - script to find longest protein coding isoform      |
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
# python ~/script/python/100b_fasta2flat.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re, classGene


### global variables
global ifile

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
    global ifile
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
            
    logfile(ifile)
            
    
def parse():
    list_obj = []
    for line in open(ifile,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            obj = classGene.GFF3(line)
        
            if obj.types() == 'mRNA':
                obj_gene = classGene.LongestCodingIsoform(line, obj)
                list_obj.append(obj_gene)
            if obj.types() == 'exon':
                obj_gene.addExon(line, obj)
            if obj.types() == 'CDS':
                obj_gene.addCDS(line, obj)
                
    hash_geneID = {}
    for obj in list_obj:
        id, parent, exon_length, cds_length = obj.getLongestIsoform()
        
        if parent in hash_geneID:
            if hash_geneID[parent][1] < cds_length:
                hash_geneID[parent][1] = cds_length
                hash_geneID[parent][2] = exon_length
                hash_geneID[parent][0] = id
            elif hash_geneID[parent][1] == cds_length:
                if hash_geneID[parent][2] < exon_length:
                    hash_geneID[parent][2] = exon_length
                    hash_geneID[parent][0] = id
        else:
            hash_geneID[parent] = [id, cds_length, exon_length]
            
    
    hash_tid = {}
    for item in hash_geneID:
        hash_tid[hash_geneID[item][0]] = ''
        
    
    ### process the GFF3 file
    print_flag = False
    for line in open(ifile,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            if obj.types() == 'gene':
                print line
            elif obj.types() == 'mRNA':
                if str(obj) in hash_tid:
                    print line
                    print_flag = True
                else:
                    print_flag = False
            elif print_flag == True:
                print line
            

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    parse()
    
    ### close the logfile
    o.close()