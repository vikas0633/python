#-----------------------------------------------------------+
#                                                           |
# 108_filterExactOverlapGFF3.py - GFF3 filtering GMAP       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# Script to filter out exact overalpping gene models        |
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
FILTER_RANGE = 50 ### remove genemodels overlapping +/- N genes

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "108_filterExactOverlapGFF3.py"+'\n')
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
            python 108_filterExactOverlapGFF3.py -i <ifile>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 108_filterExactOverlapGFF3.py -i <ifile>
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
                     

def process_file(infile):
    flag = True
    gene_coords = {}
    for line in open(infile, 'r'):
        if (len(line) > 5) or (not line.startswith('#')):
            obj = GFF3(line)
            key = obj.seqids(), int(obj.starts()), int(obj.ends())
            if obj.types() == "gene":
                if key not in gene_coords:
                    for i in range( int(obj.starts())-FILTER_RANGE, int(obj.starts())+FILTER_RANGE):
                        for j in range( int(obj.ends())-FILTER_RANGE, int(obj.ends())+FILTER_RANGE):
                            key = obj.seqids(), i, j
                            gene_coords[key] = ''
                    flag = True
                else:
                    flag = False
            if flag == True:
                print line

if __name__ == "__main__":
    
    infile = options(sys.argv[1:])
    
    process_file(infile)
    
    ### close the logfile
    o.close()