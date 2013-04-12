# 21h_fasta_AA-BP_fraction.py
# script simply takes fasta file as an input and prints
# nucleotides/AAs counts/fractions 



import os,sys, getopt
### function to calculate length of each seq in fasta file


def options(argv):
    infile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
       print 'python 21h_fasta_AA-BP_fraction.py -i <inputfile>'
       sys.exit(2)
    for opt, arg in opts:
       if opt == '-h':
          print 'python 21h_fasta_AA-BP_fraction.py -i <inputfile>'
          sys.exit()
       elif opt in ("-i", "--ifile"):
          infile = arg
    return infile
    
    
def cal_fract(infile):
	count = {}
	for line in open(infile,'r'):
		line = line.strip()
		if line[0] != '>':
			for char in line:
				if char in count:
					count[char] += 1
				else:
					count[char] = 1
	for key in sorted(count):
		print key, count[key]

if __name__ == "__main__":
    
    infile = options(sys.argv[1:])    
    
    cal_fract(infile)
    