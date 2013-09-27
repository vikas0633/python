#-----------------------------------------------------------+
#                                                           |
# 21at_FindLongestProtein.py - Finds Longest Protein        |
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
# python ~/script/python/21at_FindLongestProtein.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "21at_FindLongestProtein.py"+'\n')
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
            
    logfile(infile)
            
    return infile
    
                                
def longestIsoform(file):
    
    seqs = {}
    seqs_length = {}
    comp_header = {}
    for line in open(file):
        line = line.strip()
        
        if line.startswith('>'):
            header = line.split('.')[0]
            header_line = line
            
            
        else:
            if header in seqs_length:
                if seqs_length[header] < len(line):
                    seqs[header] = line
                    seqs_length[header] = len(line)
                    comp_header[header] = header_line
            else:
                seqs[header] = line
                seqs_length[header] = len(line)
                comp_header[header] = header_line
            
    for header in seqs:
        #print comp_header[header]
        Lj_id = comp_header[header].split(',')[0]
        try:
            Anno = ' '.join(comp_header[header].split('|')[2].split(' ')[1:]) + ' ' +comp_header[header].split('|')[0].split(',')[-1]+'|' +comp_header[header].split('|')[1]+'|'+comp_header[header].split('|')[2].split(' ')[0]
            print Lj_id +' '+ Anno
            print seqs[header]
        except:
            print comp_header[header].replace(',', ' ')
            print seqs[header]
        
            
if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    longestIsoform(file)
    
    ### close the logfile
    o.close()