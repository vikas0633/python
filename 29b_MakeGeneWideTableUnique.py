
'''
Created on June 06, 2013

@author: vgupta
'''

### script was made for making genewide table of snpEff
# 29a_MakeGeneWideTable.py - /Users/vikas0633/Desktop/script/python/ - script to put the snpEff data togehter

# Usage 
# 

import os,sys,getopt, re
import sys
import os.path
import glob

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')

### Global variables
files_prefix=[ \
"col_5", \
"col_6" ,\
"col_7", \
"col_8", \
"col_9", \
"col_10", \
"col_11", \
"col_12", \
"col_13", \
"col_14", \
"col_15", \
"col_16", \
"col_17", \
"col_18", \
"col_19", \
"col_20", \
"col_21", \
"col_22", \
"col_23", \
"col_24", \
"col_25", \
"col_26", \
"col_27", \
"col_28", \
"col_29", \
"col_30", \
"col_31", \
"col_32", \
"col_33", \
"col_34" ]

EFFECTS = ["NON_SYNONYMOUS_CODING", \
           "NON_SYNONYMOUS_START", \
           "SPLICE_SITE_ACCEPTOR", \
           "SPLICE_SITE_DONOR", \
           "START_GAINED", \
           "START_LOST", \
           "STOP_GAINED", \
           "STOP_LOST", \
           "SYNONYMOUS_CODING", \
           "SYNONYMOUS_STOP", \
           "UTR_3_PRIME", \
           "UTR_5_PRIME"]

EFFECTS_ABBREVATIONS = { "NON_SYNONYMOUS_CODING":"NSC", \
           "NON_SYNONYMOUS_START":"NSS", \
           "SPLICE_SITE_ACCEPTOR":"SSA", \
           "SPLICE_SITE_DONOR":"SSD", \
           "START_GAINED":"STG", \
           "START_LOST":"STL", \
           "STOP_GAINED":"SPG", \
           "STOP_LOST":"SPL", \
           "SYNONYMOUS_CODING":"SYC", \
           "SYNONYMOUS_STOP":"SYS", \
           "UTR_3_PRIME":"U3P", \
           "UTR_5_PRIME":"U5P"}

COLNAME="mg004	mg010	mg012	mg019	mg023	mg036	mg049	mg051	mg062	mg072	mg073	mg077	mg080	mg082	mg083	mg086	mg089	mg093	mg095	mg097	mg101	mg107	mg109	mg113	mg118	mg122	mg123	mg128	Gifu	Burtii".split()

HEADER = '#'+'Lj30_ID'+ '\t'+ 'chr'+ '\t'+ 'start'+ '\t'+ 'end'+ '\t'+ 'old_IDs'

for col in COLNAME:
    for eff in EFFECTS:
        HEADER += '\t' + col + '_' + eff
        
HEADER += '\t' + 'details'

### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    

### main argument to 
def options(argv):
    directory = ''
    suff= ''
    gene_ID = ''
    try:
        opts, args = getopt.getopt(argv,"hd:s:g:",["dir=","suff=","gene_ID="])
    except getopt.GetoptError:
        print '''
            python 29b_MakeGeneWideTableUnique.py \
                                            -d <dir> \
                                            -s <suff> \
                                            -g <gene_ID>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
            python 29b_MakeGeneWideTableUnique.py \
                                            -d <dir> \
                                            -s <suff> \
                                            -g <gene_ID>
                '''
            sys.exit()
        elif opt in ("-d", "--dir"):
            directory = arg
        elif opt in ("-s", "--suff"):
            suff = arg
        elif opt in ("-g", "--gene_ID"):
            gene_ID = arg

            
    logfile(directory)
            
    return directory, suff, gene_ID
    
### process file one by one
def process(f,snp_data):
    name = COLNAME[files_prefix.index(f.replace(suff,'')[:-1])]
    for line in open(f,'r'):
        if len(line) > 0:
            if not line.startswith('#'):
                token = line.split('\t')
                effect = token[15]
                transcript_id = token[12].strip()
                if token[0].strip() == '':
                    chro = 0
                pos = str(chro)+'_'+token[1]   
                if effect in EFFECTS:
                    if not transcript_id in snp_data:
                        snp_data[transcript_id] = {}
                    
                    if effect in snp_data[transcript_id]:
                        if pos not in snp_data[transcript_id][effect]:
                            snp_data[transcript_id][effect].append(pos)
                    else:
                        snp_data[transcript_id][effect] = [pos]
                        
                    ### adding details
                    if "details" in snp_data[transcript_id]:
                        snp_data[transcript_id]["details"] += name+'-'+str(chro) +'-'+token[1]+'-'+EFFECTS_ABBREVATIONS[effect]+'-'+token[16]+'-'+token[17]+';' 
                    else:
                        snp_data[transcript_id]["details"] = name +'-'+str(chro) +'-'+token[1]+'-'+EFFECTS_ABBREVATIONS[effect]+'-'+token[16]+'-'+token[17]+';'
    return snp_data

### Loop over the files in a folder
def LoopFiles(direcotry, suff):
    snf_data = {}
    for root, dirs, files in os.walk(direcotry):
        files = sorted(files)
        for f in files_prefix:
            f = f + "." + suff
            snf_data = process(f,snf_data)
    return snf_data   
    
### process the gene_ID information
def GeneID(gene_ID):
    gene_id_hash = {}
    for line in open(gene_ID,'r'):
        if len(line) > 1:
            if line[0] != '#':
                token = line.strip().split('\t')
                gene_id_hash[token[4]] = "\t".join(token[0:4])
    return gene_id_hash

def print_out(gene_id_hash, snp_data, suff):
    print HEADER
    for transcript_id in sorted(gene_id_hash):
        line = ''
        details = ''
        line += transcript_id + '\t'+ gene_id_hash[transcript_id]
        for effect in EFFECTS:
            if transcript_id not in snp_data:
                line += '\t' + str(0)
            else:
                if effect not in snp_data[transcript_id]:
                    line += '\t' + str(0)
                else:
                    line += '\t' + str(len(snp_data[transcript_id][effect]))
    
        if transcript_id in snp_data:
            if "details" in snp_data[transcript_id]:
                details += snp_data[transcript_id]["details"] + ','
        
        print line + '\t' + details




if __name__ == "__main__":
    direcotry, suff, gene_ID = options(sys.argv[1:])
    
    ### process the gene_ID information
    gene_id_hash = GeneID(gene_ID) 
    
    ### loop over the files in the dir
    snp_data = LoopFiles(direcotry, suff)
    
    
    ### print output from the snp_data using GENE_ID
    print_out(gene_id_hash, snp_data, suff)
        
    
    ### close the logfile
    o.close()
    
    

