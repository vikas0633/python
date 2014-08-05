#-----------------------------------------------------------+
#                                                           |
# 31h_add_FeatureType.py - script to modify GFF3 second column   |
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
global infile

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
            python 31h_add_FeatureType.py -i <ifile>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global infile
    infile = ''
    
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
            

def reaplce_feature():
    out = open(infile+'.biotype','w')
    for line in open(infile,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            tokens = line.split('\t')
            if re.search(r'Type2="NULL"',line):
                out.write(tokens[0]+'\t'+'processed_transcript'+'\t'+'\t'.join(tokens[2:])+';\n')
            elif re.search(r'Type2="repeat"',line):
                out.write(tokens[0]+'\t'+'processed_transcript'+'\t'+'\t'.join(tokens[2:])+';\n')
            elif re.search(r'Type2="protein_coding"',line):
                out.write(tokens[0]+'\t'+'protein_coding'+'\t'+'\t'.join(tokens[2:])+';\n')
            else:
                sys.exit('Correct type not found' + line)
    out.close()

if __name__ == "__main__":
    

    options(sys.argv[1:])
    
        
    reaplce_feature()
     
    o.close()
    