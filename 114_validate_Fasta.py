#-----------------------------------------------------------+
#                                                           |
# 114_validate_Fasta.py - Fasta validation script           |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# Short script to test fasta file                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/114_validate_Fasta.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables

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
            python 114_validate_Fasta.py -i <ifile>
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
    
                            
def duplicate(infile):
    headers = {}
    for line in open(infile):
        if line.startswith('>'):
            header = line[1:].strip().split(',')[0]
            if header in headers:
                print 'Error at line'
                print line
                sys.exit('Duplicate fasta header. Header is the first part of the line stariting with ">" ')
            headers[header] = ''

def findStops(file):

    for line in open(file,'r'):
        line = line.strip()
        if len(line) > 1:
            if not line.startswith('>'):
                if re.search('.',line):
                    print 'Error at line'
                    print line
                    sys.exit('Stop codon found')


if __name__ == "__main__":
    
    infile = options(sys.argv[1:])
    
    ### fasta file
    #check the duplicacy in the fasta file
    duplicate(infile)

    ### check if there is a stop in the sequence
    stops(infile)
    
    ### close the logfile
    o.close()