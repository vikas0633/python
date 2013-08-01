#-----------------------------------------------------------+
#                                                           |
# 21ao_keep_fasta_ifGFF3.py - script to throw out excessive |
# sequences in fasta file                                   |
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
# python ~/script/python/21ao_keep_fasta_ifGFF3.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "21ao_keep_fasta_ifGFF3.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 21ao_keep_fasta_ifGFF3.py
                    -i <ifile> fasta file
                    -g <gff3>  GFF3 file
            '''
    sys.exit(2)

### main argument to 

def options(argv):
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
            
    return infile, gff3
 
def empty_file(infile):
    if os.stat(infile).st_size==0:
        sys.exit('File is empty')
       
def LOADfasta(file):
    first_line = True
    seq = {}
    full_headers = {}
    string = ''
    for line in open(file,'r'):
            line = line.strip()
            if len(line) > 0 :			
                    if line[0] == '>':
                            if first_line == False:
                                    if string != '': 
                                            seq[header] = string
                            string = ''
                            header = line[1:].strip().split(',')[0]
                            full_headers[header] = line
                    else:
                            string += line
            first_line = False			
    if string != '': 
            seq[header] = string
    return seq, full_headers

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


def parseGFF3(gff3, seq, full_headers):
    for line in open(gff3, 'r'):
        line = line.strip()
        if len(line) > 1:
            obj = GFF3(line)
            
            if obj.types() == 'mRNA':
                if str(obj) in seq:
                    t_id = str(obj)
                    print full_headers[t_id]
                    print seq[t_id]
                
if __name__ == "__main__":
    
    infile, gff3 = options(sys.argv[1:])
    
    
    ### hash the fasta file
    seq, full_headers = LOADfasta(infile)
    
    ### parse GFF3 file
    parseGFF3(gff3, seq, full_headers)
    
    ### close the logfile
    o.close()