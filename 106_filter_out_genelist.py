#-----------------------------------------------------------+
#                                                           |
# 106_filter_out_genelist.py -  filter files                |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# This script takes a set of gene names and removes these   |
# genes from the file                                       |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/106_filter_out_genelist.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "106_filter_out_genelist.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    

### main argument to 

def options(argv):
    infile = ''
    format = ''
    try:
        opts, args = getopt.getopt(argv,"hi:f:g:",["ifile=","format=","list="])
    except getopt.GetoptError:
        print '''
                python 106_filter_out_genelist.py
                -i <inputfile>  
                -f <format>   [GFF3 or FASTA]
                -g <list> [file with genes, one id per line]
                '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 106_filter_out_genelist.py
                -i <inputfile>  
                -f <format>   [GFF3 or FASTA]
                -g <list> [file with genes, one id per line]
                '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-f", "--format"):
            format = arg
        elif opt in ("-g", "--list"):
            gene_ids = arg
            
    logfile(infile)
            
    return infile, format, gene_ids
    

def list_of_rows(infile):
    lst = []
    for line in open(infile):
        line = line.strip()
        if len(line) > 1:
            if line[0] != '#':
                lst.append(line)
               
    return lst
                
def filter_fasta(infile, keys):
    for line in open(infile):
        line = line.strip()
        if len(line) > 1:
            if line[0] != '#':
                if line[0] == '>':
                    for key in keys: 
                        if line[1:].startswith(key):
                            flag = False
                            break
                        else:
                            flag = True
                    
                if flag == True:
                   print line

def get_ID(line):
    match = re.search(r'ID=.+;',line)
    if match:
        return match.group().split(';')[0].replace('ID=','')
    
def filter_gff3(infile, keys):
    for line in open(infile):
        if len(line) > 1:
            if line[0] != '#':
                line = line.strip()
                token = line.split('\t')
                if token[2] == "gene":
                    id = get_ID(line)
                    if id in keys:
                        flag = False
                    else:
                        flag = True
            
                if flag == True:
                   print line
                
                            
                            
def filter_files(infile, format, keys):
    if format.lower() == "fasta":
        filter_fasta(infile, keys)
    elif format.lower() == "gff3":
        filter_gff3(infile, keys)
    else:
        sys.exit("unsupported format option")
        

if __name__ == "__main__":
    
    infile, format, gene_ids = options(sys.argv[1:])
  
    ### hash the gene_ids
    genes = list_of_rows(gene_ids)
    
    ### filter the files
    filter_files(infile, format, genes)

    ### close the logfile
    o.close()