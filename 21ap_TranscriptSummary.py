#-----------------------------------------------------------+
#                                                           |
# 21ap_TranscriptSummary.py - Summerizes GFF3 transcript wise|
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
import os,sys,getopt, re


### global variables

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

class mRNA:
    def __init__(self, line):
        self.id = get_ID(line)
        self.exon_count = 0
        self.exon_length = 0
        self.cds_count = 0
        self.cds_length = 0
    
    def __str__(self):
        return self.id
    
    def add_exon(self, line):
        token = line.split('\t')
        legnth = int(token[4]) - int(token[3]) + 1
        self.exon_count += 1
        self.exon_length += legnth

    def add_cds(self, line):
        token = line.split('\t')
        legnth = int(token[4]) - int(token[3]) + 1
        self.cds_count += 1
        self.cds_length += legnth
    
    def results(self):
        
        return self.exon_count, self.exon_length

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

def parse(gff3):
    mRNA_objs = []
    for line in open(gff3,'r'):
        if len(line) > 1 and not line.startswith('#'):
            line = line.strip()
            obj = GFF3(line)
            
            if obj.types() == 'mRNA':
                obj_mRNA = mRNA(line)
                mRNA_objs.append(obj_mRNA)
            
            if obj.types() == 'exon':
                obj_mRNA.add_exon(line)
            
            if obj.types() == 'CDS':
                obj_mRNA.add_cds(line)
        
    ### print header
    print ID +'\t'+ 'exon_count' +'\t'+ 'exon_length'
    for obj in mRNA_objs:
        exon_count, exon_length = obj.results()
        
        print str(obj) +'\t'+ str(exon_count) +'\t'+ str(exon_length)
        

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    parse(file)
    
    ### close the logfile
    o.close()