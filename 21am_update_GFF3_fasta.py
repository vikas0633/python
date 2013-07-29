#-----------------------------------------------------------+
#                                                           |
# 21am_update_GFF3_fasta.py - Update GFF3 and Fasta         |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# this script updates GFF3 and fasta given a different file |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/21am_update_GFF3_fasta.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


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
        opts, args = getopt.getopt(argv,"hi:u:f:",["ifile=","update_file=","format="])
    except getopt.GetoptError:
        print '''
            python 21am_update_GFF3_fasta.py \
                                                -i <ifile> \
                                                -u <update_file> \
                                                -f <format> 
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 21am_update_GFF3_fasta.py \
                                                -i <ifile> \
                                                -u <update_file> \
                                                -f <format> 
            '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-u", "--update_file"):
            update_file = arg
        elif opt in ("-f", "--format"):
            file_format = arg 
            
    logfile(infile)
            
    return infile, update_file, file_format
    

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
               
def process_gff3(infile):
    gene_ids = {}
    for line in open(infile, 'r'):
        if (len(line) > 1) & (not line.startswith('#')):
            obj = GFF3(line)
            if obj.types() == "gene":
                g_id = str(obj)
                if g_id in gene_ids:
                    gene_ids[g_id] += line
                else:
                    gene_ids[g_id] = line
            else:
                gene_ids[g_id] += line
            
    return gene_ids

def process_fasta(infile):
    fasta_id = {}
    for line in open(infile, 'r'):
        line = line.strip()
        if line.startswith('>'):
            header = line
            fasta_id[header] = ''
        else:
            fasta_id [header] += line
        
    return fasta_id

def hash_transcript(update_file, file_format):
    
    if file_format.lower() == 'gff3':
        return process_gff3(update_file)
    
    if file_format.lower() == 'fasta':
        return process_fasta(update_file)
    
def parse(infile, file_format, ids):
    if file_format.lower() == 'gff3':
        flag = True ### flag used for priting
        for line in open(infile,'r'):
            line = line.strip()
            if (len(line) > 1) & (not line.startswith('#')):
                obj = GFF3(line)
                if obj.types() == "gene":
                    if str(obj) in ids:
                        flag = False
                        print ids[str(obj)].strip()
                    else:
                        flag = True
            
            if flag == True:
                print line
                
    if file_format.lower() == 'fasta':
        flag = True
        for line in open(infile,'r'):
            line = line.strip()
            if line.startswith('>'):
                header = line.split(',')[0]
                if header in ids:
                    print line
                    print ids[header].strip()
                    flag = False
                else:
                    flag = True
            
            if flag == True:
                print line
            

if __name__ == "__main__":
    
    infile, update_file, file_format = options(sys.argv[1:])
    
    ### load the transcript IDs to update
    ids = hash_transcript(update_file, file_format)
    
    parse(infile, file_format, ids)
    
    
    ### close the logfile
    o.close()