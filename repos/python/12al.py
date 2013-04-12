### 12al.py - parse the mirdeep 2 output to the mysql supported output
#python 12al.py -i result_30_05_2012_t_13_17_20.csv -o miRNA_predictions

import sys, getopt


### main argument to 

def options(argv):
    inputfile = ''
    outputfile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
       print 'python 12al.py -i <inputfile> -o <outputfile>'
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print 'test.py -i <inputfile> -o <outputfile>'
          sys.exit()
       elif opt in ("-i", "--ifile"):
          inputfile = arg
       elif opt in ("-o", "--ofile"):
          outputfile = arg
    
    return inputfile, outputfile

def parse_file(inf,outf):
    parse = False
    
    def u2t(seq):
        seq_t=''
        for i in seq:
            if i == 'u':
                seq_t += 't'.upper()
            else:
                seq_t += i.upper()
        return seq_t
    
    ### open the output file
    o = open (outf,'w')
    
    for line in open(inf,'r'):
        line = line.strip()
        if (len(line) >0):
            if line.startswith('mature'):
                parse = False
                o.close()
            if line.startswith('provisional id'):
                parse = True
            if parse == True:
                token = line.split('\t')
                string=''
                for i in range(6):
                    o.write(token[i]+'\t')
                o.write(u2t(token[13])+'\t')
                for i in range(6,13):
                    o.write(token[i]+'\t')
                o.write(token[14]+'\t'+token[15]+'\n')
                
if __name__ == "__main__":
    
    inf,outf = options(sys.argv[1:])
    
    ### open the result output from the mirdeep2 and convert it to the mysql supported format
    # 1. col 7 must be the read (12al.py)
    # 2. read must only contain capital ATGC letters
    
    parse_file(inf,outf)
    
    
    
    