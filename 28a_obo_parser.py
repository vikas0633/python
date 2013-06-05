'''
Created on Feb 25, 2013

@author: vgupta
'''

### 28a_obo_parser.py - /Users/vikas0633/Desktop/script/python/ - script to obo file from the geneontology.org

### Usage
# python ~/script/python/28a_obo_parser.py -i 02_gene_ontology_ext.obo.txt -c name,namespace,def

import os,sys,getopt, re



### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile,cols):
    o.write("Program used: \t\t%s" % "28a_obo_parser.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
    o.write("Columns used: \t\t%s" % cols+'\n')        
    

### main argument to 

def options(argv):
    infile = ''
    cols = ''
    try:
        opts, args = getopt.getopt(argv,"hi:c:",["ifile=","cols="])
    except getopt.GetoptError:
        print '''
            python 28a_obo_parser.py 
                -i <inputfile>     ### obo file
                -c <cols>       ### columns to be included after GO id (comma seperated) i.e. name,namespace
                
            Example: python ~/script/python/28a_obo_parser.py -i 02_gene_ontology_ext.obo.txt -c name,namespace,def
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
            python 28a_obo_parser.py 
                -i <inputfile>     ### obo file
                -c <cols>       ### columns to be included after GO id (comma seperated) i.e. name,namespace
                
                Example: python ~/script/python/28a_obo_parser.py -i 02_gene_ontology_ext.obo.txt -c name,namespace,def
                '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-c", "--cols"):
            cols = arg            
    
    logfile(infile,cols)
            
    return infile,cols

def print_GO(file,cols):
    ''' (str, list -> '')
    this function process the file sections and print flat format
    '''
    
    first_line = True
    ### convert the cols to the list
    cols=cols.replace('"','').split(',')
        
    ### open an output file
    out = open(file+'.out','w')
    
    ### print header
    out.write("GO_id\t"+'\t'.join(cols)+'\n')

    
    ### open the file
    for line in open(file,'r'):
        line = line.strip()
        if line.startswith("[Term]"):                  ### check for new section 
            if first_line == False:
                out.write('\n')
            first_line = False
        elif line.startswith('id:'):
            token = line.split('id:')
            out.write(token[1].strip())
        else:
            token = line.split(':')
            if token[0] in cols:
                token = line.split(token[0]+':')
                out.write('\t'+token[1].strip())
        
    out.write('\n')
    out.close() ### close the output file

if __name__ == "__main__":
    
    file,cols = options(sys.argv[1:])
    
    ### col the print fucntion
    print_GO(file,cols) 
    
    ### close the logfile
    o.close()