class mRNA:
    def __init__(self, line, obj):
        self.mRNA_length = 0
        self.exon_length = 0
        self.cds_length = 0
        self.id = str(obj)
        self.hash_exon = {}
        self.hash_cds = {}
    
    def __str__(self):
        return self.id
        
    def AddData(self, line, obj):
        
        if obj.types() == 'mRNA':
            self.mRNA_length += int(obj.ends()) - int(obj.starts())
        
        if obj.types() == 'exon':
            self.exon_length += int(obj.ends()) - int(obj.starts())
            for i in range(int(obj.starts()), int(obj.ends()), 1):
                self.hash_exon[i] = ''
        
        if obj.types() == 'CDS':
            self.cds_length += int(obj.ends()) - int(obj.starts())
            for i in range(int(obj.starts()), int(obj.ends()), 1):
                self.hash_cds[i] = ''
        
    def GetExonicOverlap(self, hash): ### takes a hash with co-ordinates for exonic overlap
        return len(set(hash).intersection(set(self.hash_exon)))
    
    def GetCDSOverlap(self, hash): ### takes a hash with co-ordinates for exonic overlap
        return len(set(hash).intersection(set(self.hash_cds)))
        
    def GetmRNALength(self):
        return self.mRNA_length
        
    def GetExonLength(self):
        return self.exon_length
    
    def GetCDSLength(self):
        return self.cds_length