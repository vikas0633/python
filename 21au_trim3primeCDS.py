#-----------------------------------------------------------+
#                                                           |
# 21au_trim3primeCDS.py - script to trim 3' end of cds sequences |
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
# python ~/script/python/100b_fasta2flat.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables
global cds, protein, validate

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
            python 21au_trim3primeCDS.py -c <cds >
                                         -p <protein>
                                         -v <validate>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global cds, protein, validate
    cds = ''
    protein = ''
    validate = False
    try:
        opts, args = getopt.getopt(argv,"hc:p:v",["cds=","protein=","validate"])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-c", "--cds"):
            cds = arg
        elif opt in ("-p", "--protein"):
            protein = arg
        elif opt in ("-v", "--validate"):
            validate = True
            
    logfile(cds)
            

def hash_fasta(fasta):
    hash_headers = {}
    hash_proteins = {}
    for line in open(fasta, 'r'):
        if len(line) > 1 and not line.startswith('#'):
            line = line.strip()
            if line.startswith('>'):
                header = line.split('.')[0]
                hash_headers[header] = line
                hash_proteins[header] = ''
            else:
                hash_proteins[header] += line 
    
    return hash_proteins, hash_headers

def valid(hash_proteins):
    for line in open(cds, 'r'):
        line = line.strip()
        if len(line) > 1 and not line.startswith('#'):
            if line.startswith('>'):
                header = line.split('.')[0]
                ### check if header in the protein fasta file
                if header not in hash_proteins:
                    print 'Error at line'
                    print line
                    sys.exit('Fasta header not found in protein file')



### check if the number of bases multiple of three
def validate_cdsProteinTranslation(hash_cds, hash_headers_cds, hash_proteins, hash_headers_proteins):
    for seq in hash_cds:
        print len(hash_cds[seq]), len(hash_proteins[seq])
        if len(hash_cds[seq]) != len(hash_proteins[seq]):
            print 'Error at sequences'
            print seq
            print hash_cds[seq]
            print hash_proteins[seq]
            sys.exit('Fasta header not found in protein file')
        else:
            print seq

### check if the number of bases multiple of three
def fix_cdsProteinTranslation(hash_cds, hash_headers_cds, hash_proteins, hash_headers_proteins):
    c = open(str(now.strftime("%Y-%m-%d_%H%M."))+cds,'w')
    p = open(str(now.strftime("%Y-%m-%d_%H%M."))+protein,'w')
    for seq in hash_cds:
        c.write(hash_headers_cds[seq]+'\n')
        c.write(hash_cds[seq][:3*len(hash_proteins[seq])]+'\n')
        p.write(hash_headers_proteins[seq]+'\n')
        p.write(hash_proteins[seq]+'\n')
    c.close()
    p.close()
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### hash protein file
    hash_proteins, hash_headers_proteins = hash_fasta(protein)
    #### hash cds sequences
    hash_cds, hash_headers_cds = hash_fasta(cds)
    
    if validate == True:
        ### validate the cds file
        valid(hash_proteins)
        ### check if the number of bases multiple of three
        validate_cdsProteinTranslation(hash_cds, hash_headers_cds, hash_proteins, hash_headers_proteins)
    
    else:
        fix_cdsProteinTranslation(hash_cds, hash_headers_cds, hash_proteins, hash_headers_proteins)
        
    
    
    ### close the logfile
    o.close()