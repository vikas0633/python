#-----------------------------------------------------------+
#                                                           |
# classGeneStructure.py                                     |
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

class GeneStructure:
    
    def __init__(self, obj):
        
        ### store gene information
        self.id = str(obj)
        self.type = obj.sources()
        self.mRNAid = []
        self.mRNACoords = [] 
        self.exonid = []
        self.exonCoords = []
        self.cdsid = []
        self.cdsCoords = []
        
        ### co-ordinate hashes
        self.exons = {}
        self.cds = {}
        self.utr = {}
        self.intron = {}
        self.types = {}
    
    def __str__(self):
        return self.id
    
    def addmRNA(self, obj):
        self.mRNAid.append(str(obj))
        self.mRNACoords.append(int(obj.starts()))
        self.mRNACoords.append(int(obj.ends()))
        
    def addexon(self, obj):
        self.exonid.append(str(obj))
        self.exonCoords.append(int(obj.starts()))
        self.exonCoords.append(int(obj.ends()))
        
    def addcds(self, obj):
        self.cdsid.append(str(obj))
        self.cdsCoords.append(int(obj.starts()))
        self.cdsCoords.append(int(obj.ends()))
        
    def GetExons(self):
        self.exonCoordsSorted = sorted(self.exonCoords)
        for i in range(0, len(self.exonCoords), 2):
            for j in range(self.exonCoordsSorted[i], self.exonCoordsSorted[i+1] + 1):
                self.exons[j] = ''
        return self.exons
    
    def GetCDS(self):
        self.cdsCoordsSorted = sorted(self.cdsCoords)
        for i in range(0, len(self.cdsCoords), 2):
            for j in range(self.cdsCoordsSorted[i], self.cdsCoordsSorted[i+1] + 1):
                self.cds[j] = ''
        return self.cds
    
    def GetUTR(self):
        if len(self.cds) > 1:
            for i in self.exons:
                if i not in self.cds:
                    self.utr[i] = ''
        return self.utr
    
    def GetIntron(self):
        for i in range(0, len(self.exonCoordsSorted)-2, 2):
            for j in range(self.exonCoordsSorted[i+1] + 1 , self.exonCoordsSorted[i+2]):
                self.intron[j] = ''
        return self.intron
        
    def GetType(self):
        self.mRNACoordsSorted = sorted(self.mRNACoords) 
        for i in range(self.mRNACoordsSorted[0], self.mRNACoordsSorted[1]+1):
            self.types[i] = self.type
        return self.types