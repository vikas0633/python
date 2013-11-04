

import re

### split line
def split_line(line):
    return line.strip().split('\t')

class VCF:
    def __init__(self, line):
        
        tokens = split_line(line)
        self.CHROM = tokens[0]
        self.POS = tokens[1]
        self.ID = tokens[2]
        self.REF = tokens[3]
        self.ALT = tokens[4]
        self.QUAL = tokens[5]
        self.FILTER = tokens[6]
        self.INFO = tokens[7]
        self.FORMAT = tokens[8]
        
        self.GENOTYPE = []
        
        for i in tokens[9:]:
            self.GENOTYPE.append(i)
        
        
    def __str__(self):
        return self.CHROM+'\t'+self.POS
    
    def chroms(self):
        return self.CHROM
    
    def poss(self):
        return self.POS
    
    def ids(self):
        return self.ID
    
    def refs(self):
        return self.REF
    
    def alts(self):
        return self.ALT
    
    def quals(self):
        return self.QUAL
    
    def filters(self):
        return self.FILTER
    
    def infos(self):
        return self.INFO
    
    def formats(self):
        return self.FORMAT
    
    def genotypes(self):
        return self.GENOTYPE
    
    def depth(self):
        match = re.search(r'DP=.+;',self.INFO)
        if match:
            return match.group().split(';')[0].replace('DP=','')
        else:
            return 0
        
    def genotype(self, i):
        if len(self.GENOTYPE[i].split(':')) > 1:
            return self.GENOTYPE[i].split(':')[0]
        else:
            return 'NONE'
    
    def genotypeDepth(self, i):
        if len(self.GENOTYPE[i].split(':')) > 1:
            if len(self.GENOTYPE[i].split(':')) > self.FORMAT.split(':').index('DP') and self.genotype(i) != './.': 
                return self.GENOTYPE[i].split(':')[self.FORMAT.split(':').index('DP')]
            else:
                return 'NONE'
        else:
            return 0
    
    def genotypeQual(self, i):
        if len(self.GENOTYPE[i].split(':')) > 1:
            if len(self.GENOTYPE[i].split(':')) > self.FORMAT.split(':').index('GQ') and self.genotype(i) != './.':
                return self.GENOTYPE[i].split(':')[self.FORMAT.split(':').index('GQ')]
            else:
                return 0
        else:
            return 0
    
    def genotypeDepthSUM(self):
        geno_sum = 0
        for i in self.GENOTYPE:
            if len(i.split(':')) > 1:
                if len(self.GENOTYPE[i].split(':')) > self.FORMAT.split(':').index('DP'):
                    geno_sum += int(i.split(':')[self.FORMAT.split(':').index('DP')])
        return geno_sum
    
    def genotypeCalls(self):
        geno_call = 0
        for i in self.GENOTYPE:
            if len(i.split(':')) > 1:
                geno_call += 1
        return geno_call
    
    def genotypeCallsHete(self):
        geno_call_hete = 0
        for i in self.GENOTYPE:
            if len(i.split(':')) > 1:
                if i.split(':')[0] == '0/1' or i.split(':')[0] == '1/0':
                    geno_call_hete += 1
        return geno_call_hete
    
    def genotypeCallsHomo(self):
        geno_call_homo = 0
        for i in self.GENOTYPE:
            if len(i.split(':')) > 1:
                if i.split(':')[0] == '0/0' or i.split(':')[0] == '1/1':
                    geno_call_homo += 1
        return geno_call_homo
        
    def InbreedingCoeffs(self):
        match = re.search(r'InbreedingCoeff=.+;',self.INFO)
        if match:
            return match.group().split(';')[0].replace('InbreedingCoeff=','')
        else:
            return 0
    
    def HaplotypeScores(self):
        match = re.search(r'HaplotypeScore=.+;',self.INFO)
        if match:
            return match.group().split(';')[0].replace('HaplotypeScore=','')
        else:
            return 0
        
    def variants(self):
        if self.ALT == '.':
            return 0
        else:
            return 1
    
    def TwoPQ(self):
        zeros = 0
        ones = 0
        for i in range(len(self.GENOTYPE)):
            if self.genotype(i) == '0/0':
                zeros += 2
            elif self.genotype(i) == '0/1':
                zeros += 1
                ones += 1
            elif self.genotype(i) == '1/1':
                ones += 2
        return (zeros*ones)/float(len(self.GENOTYPE)**2)
        
        