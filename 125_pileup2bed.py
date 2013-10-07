#-----------------------------------------------------------+
#                                                           |
# 125_pileup2bed.py - Script to convert from Pileup to Bed format       |
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
            

def pileup2bed():
    last_chr = ''
    last_start = 0
    last_end = 0
    first_line = True
    name = '.'
    score = '0'
    for line in open(ifile, 'r'):
        if not line.startswith('#'):
            line = line.strip()
            token = line.split('\t')
            chro = token[0]
            pos = int(token[1])
            if first_line == False:
                if last_chr != chro and last_chr !='':
                    print str(last_chr) + '\t' + str(last_start) + '\t' + str(last_end+1) + '\t'+ name + '\t' + score
                    last_chr = chro
                    last_start = pos
                elif pos > last_end + 1 and last_chr !='':
                    print str(last_chr) + '\t' + str(last_start) + '\t' + str(last_end+1) +'\t'+ name + '\t' + score
                    last_start = pos
                if last_chr == '':
                    last_start = pos
                last_chr = chro
                last_end = pos
                
            first_line = False
    print str(last_chr) + '\t' + str(last_start) + '\t' + str(last_end + 1) +'\t'+ name + '\t' + score     
        

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    pileup2bed()
    
    ### close the logfile
    o.close()