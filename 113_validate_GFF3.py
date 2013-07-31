#-----------------------------------------------------------+
#                                                           |
# 113_validate_GFF3.py - GFF3 validation script             |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 29/07/2013                                       |
# UPDATED: 29/07/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# Script was to created to check the GFF3 file's integraty  |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/113_validate_GFF3.py -i $GFF3


### import modules
import os,sys,getopt, re, numpy


### global variables
COLS = 9

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "113_validate_GFF3.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 113_validate_GFF3.py -i <ifile> -v <verbose>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    gff3 = ''
    verbose = False
    try:
        opts, args = getopt.getopt(argv,"hi:v",["ifile=","verbose="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            gff3 = arg
        elif opt in ("-v", "--verbose"):
            verbose = True
            
    logfile(gff3)
            
    return gff3, verbose
    
                            
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
    data = {}
    gene_names = {}
    mRNA_names = {}
    for line in open(gff3,'r'):
        
        if len(line)>1:
            if not (line.startswith('#')):
                line = line.strip()
                obj = GFF3(line)
                
                ### hash gene names
                if obj.types() == "gene":
                    data[str(obj)] = line
                    if str(obj) in gene_names:
                        print 'Error at line'
                        print line
                        sys.exit('Duplicate ID'+get_ID(line))
                    else:
                        gene_names[str(obj)] = ''

                ### check orphans and hash mRNA names
                elif obj.types() == "mRNA":
                    p_ID = get_PARENT(line)
                    if p_ID not in gene_names:
                        print 'Error at line'
                        print line
                        sys.exit('Parent Missing for'+get_ID(line))
                    data[get_PARENT(line)] = line
                    if str(obj) in mRNA_names:
                        print 'Error at line'
                        print line
                        sys.exit('Duplicate ID '+get_ID(line))
                    else:
                        mRNA_names[str(obj)] = ''
                    
                else:
                    p_ID = get_PARENT(line)
                    if p_ID not in mRNA_names:
                        print 'Error at line'
                        print line
                        sys.exit('Parent Missing for '+get_ID(line))
                    
class GenomicParameters:
    
    def __init__(self, elements):
        self.key = elements.keys()[0]
        self.elements = elements
        
    def __str__(self):
        return self.key
        
    def counts(self):
        return len(self.elements[self.key])
    
    def total_length(self):
        return sum(self.elements[self.key])
    
    def avg_length(self):
        return round(numpy.mean(self.elements[self.key]),2)
    
    def max_length(self):
        return max(self.elements[self.key])
    
    def min_length(self):
        return min(self.elements[self.key])

def get_statistics(infile):
    
    print 'FeatureType'+'\t'+'Count'+'\t'+'TotalLength'+'\t'+'AverageLength'+'\t'+'MaxLength'+'\t'+'MinLength'
    gene_ids = {}
    features = {}
    elements = {'gene':[]}
    mRNA = {}
    for line in open(infile, 'r'):
        if (len(line) > 1) & (not line.startswith('#')):
            obj = GFF3(line)
            features[obj.types()]=''
            if obj.types() == "gene":
                g_id = str(obj)
                token = g_id.split('.')
                gene_id = '.'.join(token[:-1])
                if gene_id not in gene_ids:
                    gene_ids[gene_id] = ''
                elements[obj.types()].append(int(obj.ends()) - int(obj.starts()))
            else:
                ### get the mRNA length
                if obj.types() == 'mRNA':
                    mRNA_id = str(obj)
                    mRNA[mRNA_id] = 0
                if obj.types() == 'exon':
                    mRNA[mRNA_id] += int(obj.ends()) - int(obj.starts())

                
                ### other features
                if obj.types() in elements:
                    elements[obj.types()].append(int(obj.ends()) - int(obj.starts()))
                else:
                    elements[obj.types()] = [int(obj.ends()) - int(obj.starts())]
                
                
    for features in elements:
        element = {}
        element[features] = elements[features] 
        stats = GenomicParameters(element)
        
        if features == 'gene':
            print '%10s %10s %10s %10s %10s %10s' % (str(stats) , \
            str(stats.counts()) , \
            str(stats.total_length()) , \
            str(stats.avg_length()) , \
            str(stats.max_length()) , \
            str(stats.min_length()))
            gene_length = stats.total_length()
        
        
        elif features == 'mRNA':
            print '%10s %10s %10s %10s %10s %10s' % (str(stats) , \
            str(len(mRNA.values())) , \
            str(sum(mRNA.values())) , \
            str(round(numpy.mean(mRNA.values()),2)) , \
            str(max(mRNA.values())) , \
            str(min(mRNA.values())))
        
        elif features == 'exon':
            print '%10s %10s %10s %10s %10s %10s' %  (str(stats) , \
            str(stats.counts()) , \
            str(stats.total_length()) , \
            str(stats.avg_length()) , \
            str(stats.max_length()) , \
            str(stats.min_length()))
            
            print '%10s %10s %10s %10s' % (str('Intron') , \
            str(stats.counts()-len(mRNA.values())) , \
            str(gene_length - stats.total_length()) , \
            str((gene_length - stats.total_length())/(stats.counts()-len(mRNA.values()))))
        
        else:
            print '%10s %10s %10s %10s %10s %10s' %  (str(stats) , \
            str(stats.counts()) , \
            str(stats.total_length()) , \
            str(stats.avg_length()) , \
            str(stats.max_length()) , \
            str(stats.min_length()))

if __name__ == "__main__":
    
    gff3, verbose = options(sys.argv[1:])
    
    
    ### check if file empty
    empty_file(gff3)
    
    ### check if the number columns consistent
    count_cols(gff3)
    
    ### parse GFF3 file
    parse(gff3)
    
    
    ### get the statistics
    get_statistics(gff3)
    
    ### get the statistics
    print "\nYour GFF3 file has passed Vikas Gupta's completeness check\n"
    
    ### close the logfile
    o.close()
    