#-----------------------------------------------------------+
#                                                           |
# 130_shuffle_header.py - script to shuffle header of a given file |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                        |
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
import random

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
            130_shuffle_header.py  -i <ifile>
                                    -r <range> [example: 2-, 2-4]
                                    -s <sep>    [defualt: '\t']
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global infile, r, sep
    sep='\t'
    try:
        opts, args = getopt.getopt(argv,"hi:r:s:",["ifile=","range="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-r", "--range"):
            r = arg
        elif opt in ("-s", "--sep"):
            sep = arg
    
def print_out():
    
    
    start = r.split('-')[0]
    end = r.split('-')[1]
    first_line = True
    
    for line in open(infile, 'r'):
        line = line.strip()
        if len(line) > 0:
            token = line.split(sep)
            if first_line == True:
                if start == '':
                    start=0
                if end == '':
                    end = len(token)
                
                start = int(start)
                end = int(end)
                ## for indices
                start -= 1
                end -= 1
                lin = ''
                for item in token[0:int(start)]:
                    lin += item + sep
                lst = token[int(start):int(end)]
                for i in token[int(start):int(end)+1]:
                    item = random.choice(lst)
                    lst.remove(item)
                    lin += item + sep
                for item in token[int(end)+1:]:
                    lin += item + sep
                print lin[:-1]
                first_line = False
            else:
                print line
        

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### print out
    print_out()
    
    
    ### close the logfile
    o.close()