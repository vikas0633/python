#-----------------------------------------------------------+
#                                                           |
# 118_flankingRegion.py - script to extract flanking region |
# given set of co-ordinates and range                       |
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
global infile, fasta, upstream, downstream

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
            python 100b_fasta2flat.py -i <coords>
                                        -f <fasta>
                                        -u <upstream>
                                        -d <downstream>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global infile, fasta, upstream, downstream
    infile = ''
    fasta = ''
    upstream = 25
    downstream = 25
    
    try:
        opts, args = getopt.getopt(argv,"hi:f:u:d:",["coords=","fasta=","upstream=", "downstream="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--coords"):
            infile = arg
        elif opt in ("-f", "--fasta"):
            fasta = arg
        elif opt in ("-u", "--upstream"):
            upstream = int(arg)
        elif opt in ("-d", "--downstream"):
            downstream = int(arg)
            
    logfile(infile)
            
    
def Index():
    if not os.path.isfile(fasta + '.nhr'):
        os.system('formatdb -i '+fasta+' -p F -o T')

def getFlankingRegion(chro, start, end):
    os.system('nice -n 19 fastacmd -d '+ fasta+ ' -p F -s '+chro+' -L '+str(start)+','+str(end) + ' >> ' + infile+'.flanking.fasta')

def process_file():
    os.system('touch '+infile+'.flanking.fasta')
    for line in open(infile,'r'):
        if len(line)>1 and not line.startswith('#'):
            token = line.split('\t')
            chro, pos = token[0], token[1]
            start, end = int(token[1]) - upstream, int(token[1]) + downstream
            getFlankingRegion(chro, start, end)

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    ### index fasta file
    Index()
    
    ### get the co-ordinates from the file
    process_file()
    
    

    ### close the logfile
    o.close()