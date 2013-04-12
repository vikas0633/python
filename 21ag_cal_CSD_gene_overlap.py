'''
Created on March 05, 2013

@author: vgupta
'''

# 21ag_cal_CSD_gene_overlap.py - script to calculate the CDS vs gene overlap - /Users/vikas0633/Desktop/script/python

import os,sys,getopt, re
from C_loadFasta import *


### write logfile

def logfile(infile):
    import datetime
    now = datetime.datetime.now()
    o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')
    o.write("Program used: \t\t%s" % "21ag_cal_CSD_gene_overlap.py")
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
    o.close()
            
    

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print 'python 21ag_cal_CSD_gene_overlap.py -i <inputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python 21ag_cal_CSD_gene_overlap.py -i <inputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
            
    return infile


### make a fucntion to calculate the overlap
### make sure that file is sorted according to the co-ordinates and CDS follows after the genes

def cal_overlap(file): 
    first_line = True
    o = open('temp','w')
    for line in open(file,'r'):
        line = line.strip()
        token = line.split('\t')
        if token[2] == 'gene': ##chr5    GeneMark.hmm    gene    17668036    17668216    .    .    .    ID=100131_g;Name=100131_g
    
            if first_line == False:
                if len(cds) == 0:
                    cds.append(0)
                cdsSize = max(cds) - min(cds)
                #print cdsSize, GeneSize
                o.write(str(round(100*float(cdsSize)/GeneSize,2))+'\n')
            GeneSize = int(token[4]) - int(token[3]) ### get the size of the gene
            cds = []
            first_line = False
        if token[2] == 'CDS': #chr5    TAU    CDS    17668104    17668214    .    +    0    ID=100131_g.1.CDS.1;Parent=100131_g.1;
            cds.append(int(token[3]))
            cds.append(int(token[4]))
    cdsSize = max(cds) - min(cds)
    o.write(str(100*float(cdsSize)/GeneSize)+'\n')
    o.close()        
                            
        

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    cal_overlap(file)
    
    ### make the R hist plot
    os.system("R --vanilla < ~/script/R/01_plot_hist.R")