'''
Created on Feb 25, 2013

@author: vgupta
'''

### 28b_MSU_RAP_ids.py - /Users/vikas0633/Desktop/script/python/ - script to parse file containing both MSU6 and RAP-DB ids


### check output with 
# python ~/script/python/28b_MSU_RAP_ids.py -i file

import os,sys,getopt, re


### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "28b_MSU_RAP_ids.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print '''
            python 28b_MSU_RAP_ids.py -i <inputfile>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 28b_MSU_RAP_ids.py -i <inputfile>
                '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
            
    return infile
    
                            
def parse_file(file):
    """ (str) --> None
    This function will parse the ID file
    """
    for line in open(file,'r'):
        line = line.strip()
        token = line.split('\t')
        ### loop through ids in second column and print with first columns 
        for item in token[1].split(','):
            print item+'\t'+token[0]

if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    parse_file(file)
    
    ### close the logfile
    o.close()