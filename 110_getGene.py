#-----------------------------------------------------------+
#                                                           |
# 110_getGene.py - take out the genes given a list          |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# This script hash the gene_ID and take out all the mathing |
# genes from a gFF3 file                                    |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/110_getGene.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "110_getGene.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
def help():
    print '''
            python 110_getGene.py  -g <GFF3> -l <gene_list>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hg:l:",["GFF3=", "gene_list="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-g", "--GFF3"):
            gff3 = arg
        elif opt in ("-l", "--gene_list"):
            gene_list = arg
            
    logfile(gff3)
            
    return gff3, gene_list

### split line
def split_line(line):
    return line.strip().split('\t')

### get ID
def get_ID(line):
    line = line.strip()
    match = re.search(r'ID=.+;',line)
    if match:
        return match.group().split(';')[0].replace('ID=','')
    
### http://www.sequenceontology.org/gff3.shtml
### make a class that returns columns of a GFF3 row
class GFF3:
    def __init__(self, line):
        
        tokens = split_line(line)
        
        self.seqid = tokens[0]
        self.source = tokens[1]
        self.type = tokens[2]
        self.start = tokens[3]
        self.end = tokens[4]
        self.score = tokens[5]
        self.strand = tokens[6]
        self.phase = tokens[7]
        self.attribute = tokens[8]
        
        self.id = get_ID(line)
        
    def __str__(self):
        return self.id
    
    def seqids(self):
        return self.seqid
    
    def sources(self):
        return self.source
        
    def types(self):
        return self.type
    
    def starts(self):
        return self.start
    
    def ends(self):
        return self.end
    
    def scores(self):
        return self.score
    
    def strands(self):
        return self.strand
    
    def phases(self):
        return self.phase
    
    def attributes(self):
        return self.attribute

def process_objs(obj,phase, attributes):
    print obj.seqids() + '\t' + \
    obj.sources() + '\t' + \
    obj.types() + '\t' + \
    obj.starts() + '\t' + \
    obj.ends() + '\t' + \
    obj.scores() + '\t' + \
    obj.strands() + '\t' + \
    str(phase) + '\t' + \
    attributes
    
def hash_file(infile):
    hashFile = {}
    for line in open(infile,'r'):
        line = line.strip()
        hashFile[line] = ''
    
    return hashFile

def printGFF(gff3, gene_IDs):
    
    print_flag = False
    for line in open(gff3, 'r'):
        line = line.strip()
        if ((len(line) > 1) and (not line.startswith('#'))):
            obj = GFF3(line)
            if obj.types() == "gene":
                if str(obj) in gene_IDs:
                    print_flag = True
                else:
                    print_flag = False
                
            if print_flag == True:
                print line


if __name__ == "__main__":
    
    gff3, gene_list = options(sys.argv[1:])
    
    gene_IDs = hash_file(gene_list)
    
    ### print GFF
    printGFF(gff3, gene_IDs)
    
    ### close the logfile
    o.close()