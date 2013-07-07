#-----------------------------------------------------------+
#                                                           |
# 111_blastoutput_parser.py - parse blast output            |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# 111_blastoutput_parser.py - script to parse blast output  |
# and return a table                                        |
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

FASTA_FLAG = False

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
            python 111_blastoutput_parser.py 
                                -i <ifile>   ### blast output in the standard text -m6 format 
                                -f <fasta>   ### fasta file for blast database
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global FASTA_FLAG
    infile = ''
    fasta = ''
    try:
        opts, args = getopt.getopt(argv,"hi:f:",["ifile=","fasta="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-f", "--fasta"):
            fasta = arg
            FASTA_FLAG = True
            
    logfile(infile)
            
    return infile, fasta

def split_line(line, by):
    return line.split(by)

def get_species(line):
    match = re.search(r'\[(.+?)\]',line)
    if match:
        species = match.group().replace('[','').replace(']','').split(' ')
        return ' '.join(species[0:2])

def get_desc(line):
    return line.split('|')[4].split('[')[0]


### class to parse the fasta headers
class process_fasta_header:
    
    def __init__(self, line):
        
        self.token = split_line(line, '|')
        self.gi = self.token[1]
        self.desc = get_desc(line)
        self.spe = get_species(line)
        
        
    def __str__(self):
        return self.gi
    
    def GIs(self):
        return self.gi
    
    def descs(self):
        return self.desc
    
    def spes(self):
        return self.spe
    
class PARSE_BLAST:
    
    def __init__(self, line):
        
        self.token = line.split('\t')
        self.query = self.token[0]
        self.subject = '|'.join(self.token[1].split('|')[0:2])
        self.percid = self.token[2]
        self.alignlen = self.token[3]
        self.mismatch = self.token[4]
        self.gap = self.token[5]
        self.qstart = self.token[6]
        self.qend = self.token[7]
        self.sstart = self.token[8]
        self.send = self.token[9]
        self.evalue = self.token[10]
        self.bitscore = self.token[11]
        
    def __str__(self):
        return line
    
    def QUERY(self):
        return self.query
    
    def SUBJECT(self):
        return self.subject
    
    def PERCID(self):
        return self.percid
    
    def ALIGNLEN(self):
        return self.alignlen
    
    def MISMATCH(self):
        return self.mismatch
    
    def GAP(self):
        return self.gap
    
    def QSTART(self):
        return self.qstart
    
    def QEND(self):
        return self.qend
    
    def SSTART(self):
        return self.sstart
    
    def SEND(self):
        return self.send
    
    def EVALUE(self):
        return self.evalue
    
    def BITSCORE(self):
        return self.bitscore

### print headers
def header():
    
    #http://edwards.sdsu.edu/research/index.php/ramys/238-blast-output-8
    print("query "+'\t' \
          "subject" + '\t' \
          "%id" + '\t' \
          "AlignmentLength" + '\t' \
          "mismatches" + '\t' \
          "GapOpenings" + '\t' \
          "QueryStart" + '\t' \
          "QueryEnd" + '\t' \
          "SubjectStart" + '\t' \
          "SubjectEnd" + '\t' \
          "E-value" + '\t' \
          "BitScore")

### function for hashing the subject description
def hash_blast(fasta):
    blast_hash = {}
    for line in open(fasta, 'r'):
        if line.startswith('>'):
            obj = process_fasta_header(line.strip())
            
            ### hash the data against gi
            blast_hash[obj.GIs()] = [obj.descs(), obj.spes()]
    return blast_hash


### print blast results
def print_out(obj):
    
    return obj.QUERY() + '\t' + \
    obj.SUBJECT() + '\t' + \
    obj.PERCID() + '\t' + \
    obj.ALIGNLEN() + '\t' + \
    obj.MISMATCH() + '\t' + \
    obj.GAP() + '\t' + \
    obj.QSTART() + '\t' + \
    obj.QEND() + '\t' + \
    obj.SSTART() + '\t' + \
    obj.SEND() + '\t' + \
    obj.EVALUE() + '\t' + \
    obj.BITSCORE()


### function to print out blast result
def print_blast_out(infile, blast_hash):
    
    for line in open(infile,'r'):
        
        ### parse blast result
        obj = PARSE_BLAST(line.strip())
        
        ### get the blast results
        out_line = print_out(obj)
        print out_line + '\t' + blast_hash[obj.SUBJECT().split('|')[1]]
        

if __name__ == "__main__":
    
    infile, fasta = options(sys.argv[1:])
    
    ### hash the blast database
    blast_hash = hash_blast(fasta)
    
    ### print header
    header()
    
    ### parse the blast results
    print_blast_out(infile, blast_hash)
    
    ### close the logfile
    o.close()