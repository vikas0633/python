#-----------------------------------------------------------+
#                                                           |
# 115_MapFastq.py - Script to Map Fastq files               |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/115_MapFastq.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "115_MapFastq.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 115_MapFastq.py
                                    -r <ref> [reference sequence]
                                    -f <folder> [folder1, folder2, .., folderN]
                                    -x <extension> [default: fastq]
                                    -q <fastqc> [#runs fastq rather than mapping]
                                    -p <mapper> [default: bwa]
                                    -t <threads> [default: 6, numbers of core to be used]
                                    -l <seed_length> [default: 28, seed length to be used in mapping]
                                    -m <mismatches> [default: 2, mismatches allowed in the seed]
                                    -s <single> [default: Pair End alignments]
                                    -i <index> [Option to create index for Reference Sequence]
                                    -u <uncompress> [default extention: bz2]
                                    
                                    
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
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
            
    return infile
    
        

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    
    ### close the logfile
    o.close()