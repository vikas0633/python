'''
Created on Feb 25, 2013

@author: vgupta
'''

### script was made for formatting the longest ORF to flat file

### check output with 
# 21ah_count_N_between_genes.py - script to count Ns between the genes - /Users/vikas0633/Desktop/script/python
# python ~/script/python/21ah_count_N_between_genes.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame

import os,sys,getopt, re
from C_loadFasta import *
import E_get_chr_size_gff3



### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(gff3,genome):
    o.write("Program used: \t\t%s" % "21ah_count_N_between_genes.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("GFF3 file used: \t\t%s" % gff3+'\n')
    o.write("Genome file used: \t\t%s" % genome+'\n')
    

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hg:f:",["gff3=","gneome="])
    except getopt.GetoptError:
        print '''
            python 21ah_count_N_between_genes.py 
                -g <gff3> 
                -f <genome>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 21ah_count_N_between_genes.py 
                    -g <gff3> 
                    -f <genome>
                '''
        elif opt in ("-g", "--gff3"):
            gff3 = arg
        elif opt in ("-f", "--genome"):
            genome = arg
            
    logfile(gff3,genome)
            
    return gff3,genome
    
### hash positions with the Ns
def countN(genome,chr):
    for line in open(genome,'r'):
        line = line.strip()
        if line[0] == '>':
            if line[1:].split()[0] == chr:
                print line
                flag = True                 ### raise a flag if correct chromosome is found
                pos = 0
                hash = {}
                gap_count = {}
            else:
                flag = False
        elif flag == True:
            for nuc in line:                ### loop through the line
                pos += 1
                if (nuc == 'N') or (nuc == 'n'):            
                    hash[pos] = ''
                    if pos-1 not in hash:   ### a start of new gap
                        start = pos         ### save the starting positions
                        gap_count[pos]=0
                    gap_count[start] += 1
                          
                    
    ## print chr, len(hash)                 ### it prints N counts for each chromosome
    ## print chr, len(gap_count)            ### it prints number of gaps for each chromosome
    
    ''' this section print each gap length
    for key in gap_count:
        print gap_count[key]
    ''' 
    
    return hash

### this function should count number of Ns between two genes
### please remember that GFF3 file must be sorted by co-ordinates
def processGFF3(gff3,hash,chr):
    last_mRNA = 0                           ### make a stating point for chromosome
    hash = set(hash)
    for line in open(gff3,'r'):
        line = line.strip()
        if len(line) > 1:
            if line[0] != '#':
                token = line.split('\t')
                if token[2] == "mRNA":
                    match = re.search(r'ID=.+',line)
                    if match:
                        match = match.group().split(';')[0].replace('ID=','')                   
                    if token[0] == chr:
                        temp = {}           ### create a temp to store co-oridates from end of last gene till start of this gene
                        for i in range(last_mRNA,int(token[3]),1):
                            temp[i] = ''
                        print match+'\t'+str(len(set(temp).intersection(hash)))         ### this line prints the Number of Ns between two genes
                        last_mRNA=int(token[4])

if __name__ == "__main__":
    
    gff3,genome = options(sys.argv[1:])
    
    ### get chromosomes names and sizes
    chr_size = E_get_chr_size_gff3.get_size(gff3)
    
    for chr in sorted(chr_size):        ### Process file chromosomes
        hash = countN(genome,chr)       ### hash Ns in the genome by chromosomes
        processGFF3(gff3,hash,chr)
    
    ### close the logfile
    o.close()