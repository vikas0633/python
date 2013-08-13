#-----------------------------------------------------------+
#                                                           |
# 21aq_addGeneStrand.py - Adds the strand to the gene based |
# on the mRNAs strands                                      |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/21aq_addGeneStrand.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "21aq_addGeneStrand.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 21aq_addGeneStrand.py -i <ifile>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
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
            
    return infile
    
def empty_file(infile):
    if os.stat(infile).st_size==0:
        sys.exit('File is empty')
    

def count_cols(infile):
    for line in open(infile,'r'):
        if len(line)>1:
            if not (line.startswith('#')):
                if len(line.split('\t')) != COLS:
                    print 'Error at line'
                    print line
                    sys.exit('Each row must have '+str(COLS)+' columns')
                    
                    
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

def process_objs(obj, strand ):
    return obj.seqids() + '\t' + \
    obj.sources() + '\t' + \
    obj.types() + '\t' + \
    obj.starts() + '\t' + \
    obj.ends() + '\t' + \
    obj.scores() + '\t' + \
    strand + '\t' + \
    obj.phases() + '\t' + \
    obj.attributes()                        
        
class Gene:
    
    def __init__(self, obj):
        self.obj = obj
        self.features = ''
        self.strand = '.'
        
    def __str__(self):
        return str(self.obj)
    
    def add_mRNA(self, obj, line):
        if self.features == '':
            self.features = line + '\n'
        else:
            self.features += line + '\n'
        
        ### check if the strand is correct
        if self.strand != '.':
            if self.strand != obj.strands():
                pass   
                '''
                print 'Error at line'
                print line
                sys.exit('mRNAs have different headers')
                '''
        else:
            self.strand = obj.strands()
            
    def add_feature(self, obj, line):
        self.features += line + '\n'

def check_strand(file):
    
    genes = []
    for line in open(file, 'r'):
        line = line.strip()
        obj = GFF3(line)
        
        if obj.types()=="gene":
            obj_gene = Gene(obj)
            genes.append(obj_gene)
            
        elif obj.types()=="mRNA":
            obj_gene.add_mRNA(obj, line)
        
        else:
            obj_gene.add_feature(obj, line)
        
    return genes
        
def print_out(genes):
    for gene in genes:
        print process_objs(gene.obj, gene.strand)
        print gene.features

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    ### test the strand consistancy
    genes = check_strand(file)
    
    ### print out
    print_out(genes)
    
    ### close the logfile
    o.close()