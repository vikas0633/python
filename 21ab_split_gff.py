### 21ab_split_gff.py - script to split sorted GFF file based on contig/sequence/chro name

### options:
# 		split by contig
#		split by gene eveidence
#		number of contig in each file


import commands, os, sys, getopt, re, datetime
now = datetime.datetime.now()
#o = open(str(now.strftime("%Y-%m-%d_%H%M_"))+'gene.gff','w')

### Usage: python ~/script/python/21ab_split_gff.py -n 4 -i sample.gff3

### main argument to 

def options(argv):
	inputfile = ''
	splitByContig = True
	numberOfContigs = 100
	try:
		opts, args = getopt.getopt(argv,"hi:gn:",["ifile=","gene_evidence=","NumberOfContigs="])
	except getopt.GetoptError:
		print '''
			python 21ab_split_gff.py 
			-i <inputfile> 
			-g <gene_evidence> 
			-n <NumberOfContigs>
			'''
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print '''
				python 21ab_split_gff.py 
					-i <inputfile> 
					-g <gene_evidence> 
					-n <NumberOfContigs>
				'''
			sys.exit(2)
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-g", "--gene_evidence"):
			splitByContig = False
		elif opt in ("-n", "--NumberOfContigs"):
			numberOfContigs = int(arg)
	
	return inputfile, splitByContig, numberOfContigs
    
def splitGFF():
	count = -1
	last_key = ''
	
	for line in open(file):
		if len(line) > 1:
			if line[0] != '#':
				line = line.strip()
				token = line.split()
				if splitByContig == True:
					key = token[0]
				else:
					key = token[1]
				
				### eligible line 
				if key != last_key:### check if the new contig
					count += 1
					if count%numberOfContigs == 0:
						o=open(file+'_part_'+str(count/numberOfContigs + 1),'w')
				### print the line in the file
				o.write(line+'\n')
			
				last_key = key




if __name__ == "__main__":

	file, splitByContig, numberOfContigs = options(sys.argv[1:])
	
	splitGFF()
	