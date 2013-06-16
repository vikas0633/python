#-----------------------------------------------------------+
#                                                           |
# 21al_correct_strand.py - GFF3 file strand correction      |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# This script takes strand from CDS and assigns the same to |
# mRNA, exons and UTRs, GFF3 files                          |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/21al_correct_strand.py -i infile


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
            
    

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print '''
            python 21al_correct_strand.py -i <inputfile>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 21al_correct_strand.py -i <inputfile>
                '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
            
    return infile

### split line
def split_line(line):
    return line.strip().split('\t')

### get ID
def get_ID(line):
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

def process_objs(lst,strand):
    for obj in lst:
        print obj.seqids() + '\t' + \
        obj.sources() + '\t' + \
        obj.types() + '\t' + \
        obj.starts() + '\t' + \
        obj.ends() + '\t' + \
        obj.scores() + '\t' + \
        strand + '\t' + \
        obj.phases() + '\t' + \
        obj.attributes()
        
        
def parse(file):
    first_block = True
    new_block = False
    
    for line in open(file):
        obj = GFF3(line)
        if obj.types() != "gene":
            if obj.types() == "mRNA":
                if first_block == True:
                    first_block = False
                else:
                    print gene
                    process_objs(lst, strand)
                lst = []
            if obj.types() == "CDS":
                strand = obj.strands()
                
            lst.append(obj)
            
        else:
            gene = line.strip()
    
    ### for last block
    print gene
    process_objs(lst, strand)

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    parse(file)
    ### close the logfile
    o.close()