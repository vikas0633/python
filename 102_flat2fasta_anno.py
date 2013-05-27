'''
Created on March 13, 2013

@author: vgupta
'''

### script to make fasta file from the MySQL output

### check output with 
# python ~/script/python/102_flat2fasta_anno.py -i 20130313_flat.txt

import os,sys,getopt, re
from C_loadFasta import *



### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "102_flat2fasta_anno.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    

### main argument to 

def options(argv):
    infile = ''
    seq_col = 2
    anno_col = 1
    try:
        opts, args = getopt.getopt(argv,"hi:s:a:",["ifile=","seq_col=","anno_col="])
    except getopt.GetoptError:
        print '''
            python 102_flat2fasta_anno.py \
            -i <inputfile> \
            -s <seq_col> [Column containing the sequence] \
            -a <anno_col> [Column containing the annotations, separated by commas]
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 102_flat2fasta_anno.py \
                -i <inputfile> \
                -s <seq_col> [Column containing the sequence] \
                -a <anno_col> [Column containing the annotations, separated by commas]
                '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-s", "--seq_col"):
            seq_col = arg
        elif opt in ("-a", "--anno_col"):
            anno_col = arg
    
    logfile(infile)
            
    return infile, seq_col, anno_col

def makeFasta(file,seq_col,anno_col):
    first_line = True
    for line in open(file,'r'):
        line = line.strip()
        token = line.split('\t')
        if len(line)>1:
            if line[0] != '#':
                if first_line == False:
                    header = '>'
                    ''' for spider work
                    for i in list(anno_col.split(',')):
                        if (token[1] == '') or (token[1] == 0):
                            header += token[0]+','+token[0]+','
                            break   
                        header += token[int(i)-1]+','
                        '''
                    if  len(token) == int(seq_col):   
						for i in list(anno_col.split(',')): 
							if (token[1] == '') or (token[1] == '0') or (float(token[2].strip()) == 0):
								header += token[0]+','+token[-2]+','
								break
							header += token[int(i)-1]+','

						print (header.strip())[:-1]
						print token[int(seq_col)-1]
                first_line = False       
        

if __name__ == "__main__":
    
    file, seq_col, anno_col = options(sys.argv[1:])
    
    ### make fasta file
    makeFasta(file,seq_col,anno_col)
    
    
    ### close the logfile
    o.close()