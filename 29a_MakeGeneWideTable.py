
'''
Created on June 06, 2013

@author: vgupta
'''

### script was made for making genewide table of snpEff
# 29a_MakeGeneWideTable.py - /Users/vikas0633/Desktop/script/python/ - script to put the snpEff data togehter

# Usage 
# 

import os,sys,getopt, re
import sys
import os.path
import glob

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    

### main argument to 

def options(argv):
    dir = ''
    suff= ''
    try:
        opts, args = getopt.getopt(argv,"hd:s:",["dir=","suff="])
    except getopt.GetoptError:
        print '''
            python 29a_MakeGeneWideTable.py -d <dir> -s <suff>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 29a_MakeGeneWideTable.py -d <dir> -s <suff>
                '''
            sys.exit()
        elif opt in ("-d", "--dir"):
            dir = arg
        elif opt in ("-s", "--suff"):
            suff = arg
            
    logfile(dir)
            
    return dir, suff
    
                            
### Loop over the files in a folder
def LoopFiles(dir, suff):
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(suff):
                print file

if __name__ == "__main__":
    
    dir, suff = options(sys.argv[1:])
    
    ### loop over the files in the dir
    LoopFiles(dir,suff)
    
    ### close the logfile
    o.close()
    
    

