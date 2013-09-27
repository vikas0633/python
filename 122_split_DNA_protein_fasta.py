#-----------------------------------------------------------+
#                                                           |
# 122_split_DNA_protein_fasta.py - script to split fasta file between DNA and protein sequences|
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

# Example:
# python ~/script/python/122_split_DNA_protein_fasta.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables
global ifile

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 122_split_DNA_protein_fasta.py -i <ifile>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
            
    logfile(ifile)
            

def parse_fasta():
    dna_char = 'ATGCNXatgcnx'
    first_line = True
    protein = False
    d = open(ifile+'.dna','w')
    p = open(ifile+'.prot','w')
    for line in open(ifile,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            if line.startswith('>'):
                if first_line == False:
                    if protein == True:
                        p.write(header+'\n')
                        p.write(seq+'\n')
                    else:
                        d.write(header+'\n')
                        d.write(seq+'\n')
                protein = False
                header = line
                seq = ''
                first_line = False
            else:
                seq += line
                for char in line:
                    if char not in dna_char:
                        protein = True
    
    if protein == True:
        p.write(header+'\n')
        p.write(seq+'\n')
    else:
        d.write(header+'\n')
        d.write(seq+'\n')
    
    d.close()
    p.close()

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### load fasta file
    parse_fasta()
    
    ### close the logfile
    o.close()