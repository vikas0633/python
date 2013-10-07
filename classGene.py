### Gene class 

import re

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
        
        self.line = line
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
    
    def get_parent(self):
        return get_PARENT(self.line)

def process_objs(obj):
    return obj.seqids() + '\t' + \
    obj.sources() + '\t' + \
    obj.types() + '\t' + \
    obj.starts() + '\t' + \
    obj.ends() + '\t' + \
    obj.scores() + '\t' + \
    obj.strands() + '\t' + \
    obj.phases() + '\t' + \
    obj.attributes()
    
    
class Gene:
    
    def __init__(self, obj):
        self.obj = obj
        self.features = ''
        self.strand = '.'
        
    def __str__(self):
        return str(self.obj)
    
    def add_mRNA(self, obj, line):
        if self.features == '':
            self.features = line + '\n'
        else:
            self.features += line + '\n'
        
        ### check if the strand is correct
        if self.strand != '.':
            if self.strand != obj.strands():
                pass   
                '''
                print 'Error at line'
                print line
                sys.exit('mRNAs have different headers')
                '''
        else:
            self.strand = obj.strands()
            
    def add_feature(self, obj, line):
        self.features += line + '\n'

class LongestIsoform:
    def __init__(self, obj):
        self.obj = obj
        self.gene = process_objs(obj)
        self.features = ''
        self.longest_isoform_id = ''
        self.longest_isoform_length = 0
        self.PrintFeature = False
        
    def __str__(self):
        return str(self.gene)
    
    def add_mRNA(self, line, obj):
        
        if self.longest_isoform_length < int(obj.ends()) - int(obj.starts()):
            self.longest_isoform_length = int(obj.ends()) - int(obj.starts())
            self.PrintFeature = True
        else:
            self.PrintFeature = False
            
        if self.PrintFeature == True:
            self.features = line + '\n'
            
    def add_feature(self, line):
        if self.PrintFeature == True:
            self.features += line + '\n'

def check_strand(file):
    
    genes = []
    for line in open(file, 'r'):
        line = line.strip()
        obj = GFF3(line)
        
        if obj.types()=="gene":
            obj_gene = Gene(obj)
            genes.append(obj_gene)
            
        elif obj.types()=="mRNA":
            obj_gene.add_mRNA(obj, line)
        
        else:
            obj_gene.add_feature(obj, line)
        
    return genes
        
def print_out(genes):
    for gene in genes:
        print process_objs(gene.obj, gene.strand)
        print gene.features
        
class LongestCodingIsoform:
    def __init__(self, line, obj):
        self.obj = obj
        self.id = str(obj)
        self.parent = get_PARENT(line)
        self.features = ''
        self.exon_length = 0
        self.cds_length = 0
        self.data = line + '\n'
        
    def __str__(self):
        return str(self.gene)
    
    def addExon(self,line, obj):
        self.exon_length += int(obj.ends()) - int(obj.starts()) + 1
        
    def addCDS(self,line, obj):
        self.cds_length += int(obj.ends()) - int(obj.starts()) + 1
        
    
    def getLongestIsoform(self):
        return self.id, self.parent, self.exon_length, self.cds_length
    