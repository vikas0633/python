#-----------------------------------------------------------+
#                                                           |
# 16b_blast_reciprocal.py - script to parse vcf format file |
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
            python 100b_fasta2flat.py -i <ifile>
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
            
def parseBlastOut():
    hash_blast = {}
    for line in open(ifile,'r'):
        line = line.strip()
        if len(line)>0 and not line.startswith('#'):
            tokens = line.split('\t')
            hash_blast[tokens[0]] = tokens[1]
        
    pairs_out = {}
    ### look for reprocals
    for line in open(ifile,'r'):
        line = line.strip()
        if len(line)>0 and not line.startswith('#'):
            tokens = line.split('\t')
            if tokens[1] in hash_blast:
                if hash_blast[tokens[1]] == tokens[0] and tokens[1] != tokens[0]:
                    key=tokens[0]+'-'+tokens[1]
                    if key not in pairs_out:
                        print line
                        key=tokens[0]+'-'+tokens[1]
                        pairs_out[key] = ''
                        key=tokens[1]+'-'+tokens[0]
                        pairs_out[key] = ''

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    
    parseBlastOut()
    
    ### close the logfile
    o.close()