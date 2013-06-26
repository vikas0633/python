#-----------------------------------------------------------+
#                                                           |
# 109_AddPhaseGFF3.py - Adds the GFF3 phase                    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 26/06/2013                                       |
# UPDATED: 26/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# Given CDS sequences in the GFF3 file this scripts adds the|
# phase based on the exon co-ordinates                      |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/1108_AddPhaseGFF3 -i GFF3


### import modules
import os,sys,getopt, re


### global variables
CDS = {}
mRNA_id = ''
CDS_strand = ''

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
            python 108_AddPhaseGFF3-i <ifile>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                108_AddPhaseGFF3 -i <ifile>
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

def getPhase(obj):
    global CDS
    CDS[int(obj.starts())] = obj

### returns a list containing phase values
def calcPhase(CDS, start_count):
    phase = []
    if start_count == 0:
        for start in sorted(CDS.keys()):
            end = CDS[start].ends()
            if len(phase) == 0:
                phase.append(0)
                last_phase = phase[-1]
                phase.append( ( 3 -  ((int(end) - int(start) + 1 )%3) + last_phase )%3)
                last_phase = phase[-1]
            else:
                phase.append( ( 3 -  ((int(end) - int(start) + 1 )%3) + last_phase )%3)
                last_phase = phase[-1]
        return phase
    else:
        lst = sorted(CDS.keys())
        for start in lst[::-1]:
            end = CDS[start].ends()
            if len(phase) == 0:
                phase.append(0)
                last_phase = phase[-1]
                phase.append( (((int(end) - int(start) + 1 )%3) + last_phase )%3)
                last_phase = phase[-1]
            else:
                phase.append( (((int(end) - int(start) + 1 )%3) + last_phase )%3)
                last_phase = phase[-1]
        return phase[::-1][1:]

def process_CDS(count):
    
    start_count = count
    phase = calcPhase(CDS, start_count)
    i = 0
    for start in sorted(CDS.keys()):
        if start_count == 0:
            count += 1
            attributes =  re.sub(r'ID.+?;','ID='+mRNA_id+'.CDS.'+str(count)+';',CDS[start].attributes())
            process_objs(CDS[start], phase[i], attributes)
        else:
            attributes =  re.sub(r'ID.+?;','ID='+mRNA_id+'.CDS.'+str(count)+';',CDS[start].attributes())
            process_objs(CDS[start], phase[i], attributes)
            count -= 1
        i += 1   
            
def printCDS():
    if CDS_strand == '+':
        count = 0
        process_CDS(count)
    else:
        count = len(CDS) 
        process_CDS(count)

def parseGFF3(gff3):
    global CDS, mRNA_id, CDS_strand
    CDS_flag = False
    for line in open(gff3,'r'):
        if (len(line)>2) & (not line.startswith('#')): 
            obj = GFF3(line)
            if obj.types() != 'CDS':
                if CDS_flag == True:
                    printCDS()
                    CDS = {}
                    CDS_flag = False
                print line.strip()
            
            else:
                CDS_flag = True
                CDS_strand = obj.strands()
                getPhase(obj)
            
            if obj.types() == 'mRNA':
                mRNA_id = str(obj)
                
        first_line = False
    
    printCDS()
        

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    parseGFF3(file)
    
    ### close the logfile
    o.close()