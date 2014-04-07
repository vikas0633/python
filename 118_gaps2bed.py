#-----------------------------------------------------------+
#                                                           |
# 118_gaps2bed.py - script takes a fasta file and created bed file with gap co-ordinates |
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
            python 118_gaps2bed.py
                -i <ifile> ## fasta file
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
            
def make_bed():
    '''
    # bed format
    chloro	20419	20419	L0843	13	+	
    chloro	28071	28071	L5618	21	-
    '''
    
    
    o = open(infile+'.gaps.bed','w')
    for line in open(infile, 'r'):
        line = line.strip()
        
        if len(line) >0 and not line.startswith('#'):
            if line.startswith('>'):
                chro = line[1:]
                pos = 0
                start_N = 0
                end_N = 0
                N_region = False
            else:            
                for char in line:
                    pos += 1
                    if pos%1000 == 0:
                        print 'Processing positions: ', chro, '{:10,.0f}'.format(pos)
                    if char=='N':
                        if N_region == True:
                            end_N += 1
                        else:
                            start_N = pos
                            end_N = pos
                            N_region = True
                    else:
                        if N_region == True:
                            N_region = False
                            name=chro+'_'+str(start_N)+'_'+str(end_N)
                            o.write(chro+'\t'+str(start_N)+'\t'+str(end_N)+'\t'+name+'\t'+'.'+'\t'+'+'+'\n')
    o.close()
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    make_bed()
    
    ### close the logfile
    o.close()