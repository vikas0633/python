#-----------------------------------------------------------+
#                                                           |
# 21ak_remove_redundant.py - template for python scripting  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 06/10/2013                                       |
# UPDATED: 06/10/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
#  Script to remove redundant gene models                   |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/21ak_remove_redundant.py 


import os,sys,getopt, re
import E_get_chr_size_gff3

### Global variables
IGNORE="CUFF"
MIN_OVERLAP=0.8


### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')

### import modules
import os,sys,getopt, re


### global variables


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
            python 21ak_remove_redundant.py -i <inputfile>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 21ak_remove_redundant.py -i <inputfile>
                '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
            
    return infile

def get_ID(line):
    match = re.search(r'ID=.+;',line)
    if match:
        return match.group().split(';')[0].replace('ID=','')
    
def get_NAME(line):
    match = re.search(r'Name=.+;',line)
    if match:
        return match.group().split(';')[0].replace('Name=','')
    
def hash_coords(line):
    temp = {}
    token = line.split('\t')
    for i in range(int(token[3]),int(token[4])+1):
        temp[i] = ''
    return set(temp)

### GeneInfo class
class GeneInfo:
    def __init__(self, line):
        
        self.line = line
        
        ### get the id
        self.id = get_ID(line)
        
        ### get the name
        self.name = get_NAME(line)
        
        ### hash the co-ordinates
        self.coord = hash_coords(line)
    
        
    def __str__(self, line):
        
        return self.id
    
    def ids(self):
        return self.id
    
    def names(self):
        return self.name
    
    def coords(self):
        return self.coord

def parse_geneModels(f, chro):
    lst = []
    for line in open(f,'r'):
        line = line.strip()
        token = line.split('\t')
        
        if (token[0]==chro) & (token[2]=="gene"):
            lst.append(GeneInfo(line))
            
                            
    return lst


def detect_overlap(lst):
    for item in lst:
        ### get the name of the gene
        name = item.names()
        id = item.ids()
        coord = item.coords()
        if not name.startswith(IGNORE): 
            for item2 in lst[lst.index(item):]:
                if id != item2.ids():
                    if (len(coord.intersection(item2.coords())) > MIN_OVERLAP*(len(coord))):
                        print id, name, item2.ids(), item2.names()
                


if __name__ == "__main__":
    
    f = options(sys.argv[1:])
     
    size = E_get_chr_size_gff3.get_size(f)
    
    for chro in sorted(size):
        ### process the gene model file
        lst = parse_geneModels(f, chro)
        ### find overlapping gene models
        detect_overlap(lst)
    
    ### close the logfile
    o.close()