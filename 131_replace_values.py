#-----------------------------------------------------------+
#                                                           |
# 131_replace_values.py - script to replace all the values in a row  |
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
import os,sys,getopt, re, random

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
                        131_replace_values.py  -i <ifile>
                                    -r <range> [example: 2-2, 2-4, rows to be replaced]
                                    -s <sep>    [defualt: '\t']
                                    -w <word>   [default: NA, word to be used for replacing]
                                    -c <col>    [default: 2-, replace these columns]
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global infile, r, sep, word, col
    word='NA'
    col='2-'
    sep='\t'
    try:
        opts, args = getopt.getopt(argv,"hi:r:s:w:c:",["ifile=","range=","sep=","word=","col="])
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
        elif opt in ("-w", "--word"):
            word = arg
        elif opt in ("-c", "--col"):
            col = arg    
    
                            
def print_out():
    row_count= 0
    st = int(r.split('-')[0])
    en = int(r.split('-')[1])
    for line in open(infile, 'r'):
        line = line.strip()
        row_count += 1
        if len(line) > 0:
            if int(st) <= row_count <= int(en):
                
                start = col.split('-')[0]
                end = col.split('-')[1]
                token = line.split(sep)
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
                    item = word
                    lin += item + sep
                for item in token[int(end)+1:]:
                    lin += item + sep
                print lin[:-1]
            else:
                print line
                
                
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    print_out()
    
    ### close the logfile
    o.close()