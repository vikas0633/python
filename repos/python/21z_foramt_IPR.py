### 21z_foramt_IPR.py - script takes raw output from IPRScan and make non-redundant gene_ID\annotation  - /Users/vikas0633/Desktop/script/python


import os,sys, getopt


def options(argv):
	ipr = ''
	try:
		opts, args = getopt.getopt(argv,"hi:",["ipr="])
	except getopt.GetoptError:
		print 'python 21z_foramt_IPR.py -i <ipr>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'python 21z_foramt_IPR.py -i <ipr>'
			sys.exit()
		elif opt in ("-i", "--ipr"):
			infile = arg
	return infile



### parsing the IPRScan output 
### making gene-model wise annotation list
def parse(infile):
	import re
	IPR_out = {}
	for line in open(infile,'r'):
		line = line.strip()
		token = line.split('\t')
		if token[0] not in IPR_out:
			IPR_out[token[0]] = token[5]+','+token[12]
		else: ### make sure that there are not duplicate annotations
			list = IPR_out[token[0]].split(';')
			flag = True
			for i in list:
				if (token[5]+','+token[12]) == i.strip():
					  flag = False
			if flag==True: ### check if coulmn 6 or column 13 annotations are already have been included 
					IPR_out[token[0]] += '; '+token[5]+','+token[12]
			 
			
	return IPR_out

### print file
def printout(IPR):
	for key in IPR:
		print key +'\t'+ IPR[key]


if __name__ == "__main__":
	
	infile = options(sys.argv[1:])    
	
	
	### parse the IPRScan output
	IPR = parse(infile)

	### print file
	printout(IPR)
	
