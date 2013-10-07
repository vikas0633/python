#-----------------------------------------------------------+
#                                                           |
# classMAF.py - script to parse maf format file             |
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


def Coords(line):
    token = line.split()
    chro = token[0]
    
    strand = token[3]
    size = int(token[2])
    seq = token[5]
    
    TrueSize = len(seq.replace('-',''))
    
    if token[3] == '+':
        start = int(token[1])
        end = int(token[1]) + TrueSize
    if token[3] == '-':
        start = int(token[1]) - TrueSize
        end = int(token[1])
        
        start = int(token[4]) - start
        end = int(token[4]) - end
        
    
    return chro, start, end, strand, size, seq, TrueSize
    

### import modules
import os,sys,getopt, re

class MAF:
    def __init__(self, obj_count):
        self.count = obj_count
        self.data = ''
        
    def __str__(self):
        return self.count
    
    def addData(self, line):

        line = line.strip()
        self.data += '\n'+ line
        
        match = re.search(r'# identity=.+',line)
        if match:
            self.ident = match.group().split('(')[1].split(')')[0].replace('%','')
        
        match = re.search(r'# coverage=.+',line)
        if match:
            self.cover = match.group().split('(')[1].split(')')[0].replace('%','')
            
        match = re.search(r'# continuity=.+',line)
        if match:
            self.conti = match.group().split('(')[1].split(')')[0].replace('%','')
        
        match = re.search(r'# cigar=.+',line)
        if match:
            self.cig = match.group().replace('# cigar=','')
        
        match = re.search(r'a score=.+',line)
        if match:
            self.sco = match.group().replace('a score=','')
    
    def addRefAlign(self,line):
        line = line.strip()
        self.data += '\n'+ line
        
        line = line[2:]
        
        
        
        self.ref_chro, self.ref_start, self.ref_end, self.ref_strand, self.ref_size, self.ref_seq, self.ref_TrueSize = Coords(line)
        
    def addtargetAlign(self,line):
        line = line.strip()
        self.data += '\n'+ line
        
        line = line[2:]
        self.target_chro, self.target_start, self.target_end, self.target_strand, self.target_size, self.target_seq, self.target_TrueSize  = Coords(line)
    
        ### make the alignment against lotus
        ### remove the '-' positions from the Lotus
        self.alignHash = {}
        count = 0
        for i in range(len(self.ref_seq)):
            if self.ref_seq[i] != '-':
                count += 1
                self.alignHash[count] = str(self.ref_seq[i])+str(self.target_seq[i])
    
    
    def Data(self):
        return self.data
    
    def identity(self):
        return self.ident
        
    def coverage(self):
        return self.cover
    
    def continuity(self):
        return self.conti
    
    def cigar(self):
        return self.cig

    def score(self):
        return self.sco
    
    def RefChro(self):
        return self.ref_chro
    
    def RefStart(self):
        return self.ref_start
    
    def RefEnd(self):
        return self.ref_end
    
    def RefStrand(self):
        return self.ref_strand
        
    def RefSize(self):
        return self.ref_size
    
    def RefTrueSize(self):
        return self.ref_TrueSize
    
    def TargetChro(self):
        return self.target_chro
    
    def TargetStart(self):
        return self.target_start
    
    def TargetEnd(self):
        return self.target_end
    
    def TargetStrand(self):
        return self.target_strand

    def TargetSize(self):
        return self.target_size
    
    def TargetTrueSize(self):
        return self.target_TrueSize
    
    def AlignHash(self):
        return self.alignHash
    
    