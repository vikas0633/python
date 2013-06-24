#-----------------------------------------------------------+
#                                                           |
# 107_parameterGFF3.py - GFF3 parameters                    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 22/06/2013                                       |
# UPDATED: 22/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# Script to calculate Gene, mRNA, exon, CDS count, Total    |
# length and average length                                 |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/100b_fasta2flat.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re, numpy


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
            python 100b_fasta2flat.py -i <ifile>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 100b_fasta2flat.py -i <ifile>
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

           
def process_gff3(infile):
    
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
            print str(stats) + '\t' + \
            str(stats.counts()) + '\t' + \
            str(stats.total_length()) + '\t' + \
            str(stats.avg_length()) + '\t' + \
            str(stats.max_length()) + '\t' + \
            str(stats.min_length())
            gene_length = stats.total_length()
        
        
        elif features == 'mRNA':
            print str(stats) + '\t' + \
            str(len(mRNA.values())) + '\t' + \
            str(sum(mRNA.values())) + '\t' + \
            str(round(numpy.mean(mRNA.values()),2)) + '\t' + \
            str(max(mRNA.values())) + '\t' + \
            str(min(mRNA.values()))
        
        elif features == 'exon':
            print str(stats) + '\t' + \
            str(stats.counts()) + '\t' + \
            str(stats.total_length()) + '\t' + \
            str(stats.avg_length()) + '\t' + \
            str(stats.max_length()) + '\t' + \
            str(stats.min_length())
            
            print str('Intron') + '\t' + \
            str(stats.counts()-len(mRNA.values())) + '\t' + \
            str(gene_length - stats.total_length()) + '\t' + \
            str((gene_length - stats.total_length())/(stats.counts()-len(mRNA.values())))
        
        else:
            print str(stats) + '\t' + \
            str(stats.counts()) + '\t' + \
            str(stats.total_length()) + '\t' + \
            str(stats.avg_length()) + '\t' + \
            str(stats.max_length()) + '\t' + \
            str(stats.min_length())
            
    
        

if __name__ == "__main__":
    
    infile = options(sys.argv[1:])
    
    
    ### process the file
    process_gff3(infile)
    
    ### close the logfile
    o.close()