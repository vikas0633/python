### 21u_make_gff2.py - script makes gff2 file for the TAU input, same as Stig's 26_parse.pl - /Users/vikas0633/Desktop/script/python


import os,sys, getopt

def file_empty(file):
    count = sum([1 for line in open(file)])
    if count == 0:
        sys.exit(file+' is empty')	



def options(argv):
	infile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["in=","out="])
	except getopt.GetoptError:
		print 'python 21u_make_gff2.py -i <in> -o <out>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'python 21u_make_gff2.py -i <in> -o <out>'
			sys.exit()
		elif opt in ("-i", "--in"):
			infile = arg
		elif opt in ("-o", "--out"):
			outfile = arg
	
	return infile, outfile
	


def process_file():
	o = open('temp.gff','w')
	for line in open(infile,'r'):
	        line = line.strip()
		token = line.split('\t')
		if (token[2] == "exon"):
		        lin = token[0]+'\t'+token[1]+'\t'+token[2]+'\t'+token[3]+'\t'+token[4]+'\t'+token[5]+'\t'+token[6]+'\t'+token[7]+'\t'+line.split('=')[2].split(';')[0]
		        o.write(lin+'\n')
	o.close()
if __name__ == "__main__":

	infile, outfile = options(sys.argv[1:]) 
	
	
	#file_empty(infile)
	#file_empty(gff)
	
	
	### process gene model file
	process_file()