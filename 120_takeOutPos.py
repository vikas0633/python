#-----------------------------------------------------------+
#                                                           |
# 120_takeOutPos.py - script to take out matching coords    |
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
global ifile, coords

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
            python 100b_fasta2flat.py -i <ifile> -c <coords>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile, coords
    ifile = ''
    coords = ''
    try:
        opts, args = getopt.getopt(argv,"hi:c:",["ifile=","--coords="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-c", "--coords"):
            coords = arg
            
    logfile(ifile)
            
    
def hashCoords():
    hash_coords ={}
    for line in open(coords,'r'):
        if len(line) > 1 and not line.startswith('#'):
            tokens = line.split('\t')
            hash_coords[tokens[0],tokens[1]] = ''
    return hash_coords    

def parse(hash_coords):
    for line in open(ifile, 'r'):
        if len(line) > 1 and not line.startswith('##'):
            line = line.strip()
            tokens = line.split('\t')
            key = tokens[0],tokens[1]
            
            if key in hash_coords:
                print line
            if line.startswith('#'):
                print line

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    hash_coords = hashCoords()
    
    parse(hash_coords)
    
    ### close the logfile
    o.close()