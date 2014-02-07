#-----------------------------------------------------------+
#                                                           |
# 27_TranscriptsOnScaffold.py -                             |
# Script to extract all the transcripts on given scaffolds  |
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
# python ~/script/python/27_TranscriptsOnScaffold.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re
import classGene

### global variables
global infile, gff3, FeatureType

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
            python 100b_fasta2flat.py -i <ifile> -g <gff3> -t <type>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global infile, gff3, FeatureType
    infile = ''
    gff3 = ''
    FeatureType = 'mRNA'
    
    try:
        opts, args = getopt.getopt(argv,"hi:g:t:",["ifile=","gff3=","type="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-g", "--gff3"):
            gff3 = arg
        elif opt in ("-t", "--type"):
            FeatureType = arg        
        
        
    logfile(infile)
    
def HashCol():
    Hash = {}
    for line in open(infile, 'r'):
        line = line.strip()
        token = line.split()
        Hash[token[0]] = ''
    
    return Hash

def printOut(Hash):
    
    
    for line in open(gff3, 'r'):
        if len(line) > 0 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            if obj.types() == FeatureType:
                if obj.seqids() in Hash:
                    print obj

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    
    ### hash the scaffolds
    Hash = HashCol()
    
    ### print out the feature ids
    printOut(Hash)
    
    ### close the logfile
    o.close()