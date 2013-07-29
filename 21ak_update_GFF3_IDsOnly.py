#-----------------------------------------------------------+
#                                                           |
# 21ak_update_GFF3_IDsOnly.py this script take a two column |
# id and replaces these in the GFF3 file                    |
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
            python 100b_fasta2flat.py -i <ifile> -l <gene_list>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:l:",["ifile=","list="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-l", "--list"):
            list_file = arg
        
    logfile(infile)
            
    return infile, list_file
    
def hash_file(infile):
    hash_ID = {}
    for line in open(infile,'r'):
        line = line.strip()
        token = line.split('\t')
        hash_ID[token[0]] = token[1] ### hash trasctipt ID
        g_id = '.'.join(token[0].split('.')[:-1])
        g_id_new = '.'.join(token[1].split('.')[:-1])
        hash_ID[g_id] = g_id_new ### hash gene_ID
    return hash_ID

### split line
def split_line(line):
    return line.strip().split('\t')

### get ID
def get_NAME(line):
    line = line.strip()
    match = re.search(r'Name=.+',line)
    if match:
        return match.group().split(';')[0].replace('Name=','')
    
### get Name
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


def update_GFF3(gff3, hash_ID):
    
    for line in open(gff3, 'r'):
        line = line.strip()
        if ((len(line) > 1) and (not line.startswith('#'))):
            obj = GFF3(line)
            
            if obj.types() == 'gene':
                name = get_NAME(line)
                g_id = get_ID(line)
                
            if name in hash_ID:    
                print line.replace(g_id, hash_ID[name])
            else:
                print line.replace(g_id, name)
                
if __name__ == "__main__":
    
    gff3, list_file = options(sys.argv[1:])
    
    ### hash file
    hash_ID = hash_file(list_file)
    
    ### update GFF3
    update_GFF3(gff3, hash_ID)
    
    ### close the logfile
    o.close()