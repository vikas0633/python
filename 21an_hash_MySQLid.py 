
#-----------------------------------------------------------+
#                                                           |
# 21an_hash_MySQLid.py - this script makes a 2 column table |
#    one with Id and another with yes/no                    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# this script makes a 2 column table one with Id and another|
# with yes/no                                               |
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
        elif opt in ("-d", "--data"):
            data = arg
        
            
    logfile(infile)
            
    return infile, data
    
def hash_file(file):
    hash_data ={}
    for line in open(file, 'r'):
        line = line.strip()
        hash_data[line] = ''
                            
    return hash_data

def printOut(hash_data, data):
    
    for line in open(data,'r'):
        line = line.strip()
        
        if line in hash_data:
            print line + '\t' + 'yes'
        else:
            print line + '\t' + 'no'

if __name__ == "__main__":
    
    Ids, data  = options(sys.argv[1:])
    
    hash_data = hash_file(Ids)
    
    printOut(hash_data, data)
    
    ### close the logfile
    o.close()