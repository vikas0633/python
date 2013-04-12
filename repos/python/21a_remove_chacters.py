
#21a_remove_chacters.py - this script removes the any other character than ATGCN

import os,sys, getopt


### main argument to 

def options(argv):
    inputfile = ''
    outputfile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
    except getopt.GetoptError:
       print 'python 21a_remove_chacters.py -i <inputfile> -o <outputfile>'
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
	def replace_by_N(string):
		new_string=''
		for i in string:
			if i not in 'ATGCNatgcn':
				new_string += 'N'
			else:
				new_string += i
		return new_string
	### open the output file
	o = open (outf,'w')
	for line in open(inf,'r'):
		line = line.strip()
		if line.startswith('>'):
			o.write(line+'\n')
		else:
			o.write(replace_by_N(line)+'\n')  
	o.close()
                 
if __name__ == "__main__":
    
    inf,outf = options(sys.argv[1:])
    if outf == '':
    	outf = inf+'_corrected'
    
    ### parse file with the sequence
    parse_file(inf,outf)
