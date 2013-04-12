### script 12ai.py
### this script adds 3' demo adapter to a sequences with a given degree of base calling error

### Usages: nice -n 19 python 12ai.py <infile> <adaptor_sequence> <base_calling_error(0,1)>> <outfile>

### base quality of adapter is the same as the last base at the 3' end

import sys


### function for generating adapter and writing reads to the file
def write_reads(infile, adapter,bce=0.05):
    import random
    
    ### open and read the infile
    lc = 0
    for line in open(infile,'r'):
        
        ### randomly select base to be changed
        c=random.randrange(4)
        nuc='ATGC'
        
        ### make new adapter
        new_adapter=''
        for i in adapter:
            if random.random() < bce:
                c=random.randrange(4)
                new_adapter += nuc[c]
            else:
                new_adapter += i
        
        
        line = line.strip()
        lc += 1
        
        if lc%4 == 1:
            print line
        elif lc%4 == 2:
            print line+new_adapter
        elif lc%4 == 3:
            print line
        elif lc%4 == 0:
            fake_qual = [line[-1] for i in adapter]
            fake_qual=''.join(fake_qual)
            print line+fake_qual
    




if __name__ == "__main__":
    
    ### read adapter sequence
    adapter=sys.argv[2]
    
    ### read the input file
    infile=sys.argv[1]
    
    ### base calling error
    bce=float(sys.argv[3])
    
    ### call the function to write the adapter added files
    write_reads(infile,adapter,bce)