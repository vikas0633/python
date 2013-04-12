# 21ad_makebed.py - script to make bed format file from the given column names - /Users/vikas0633/Desktop/script/python

import os,sys,getopt, re

### Usage: python ~/Desktop/script/python/21ad_makebed.py -i dk012357_insertion_table_sorted.txt -c 3,4,4,9,8,5 -r '' > Lore1.bed

### main argument to 

def options(argv):
	file1 = ''
	col1 = ''
	sep='\t'
	first_line = False
	header = 'Chrom\tChromStart\tChromEnd\tName\tAbundance\tstrand'
	try:
		opts, args = getopt.getopt(argv,"hi:c:s:r:",["file=","col="'separatedBy=',"header="])
	except getopt.GetoptError:
		print '''
			python 21ad_makebed.py
				-i <file> # file containing co-ordinate information
				-c <col> # multiple columns separated by commas in same order as of bed format
				-s <separatedBy>
				-r <header> # header to be used in the bed formatted file
				-f <first_line>
				-h <help> 
			'''
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print '''
				python 21ad_makebed.py
					-i <file> # file containing co-ordinate information
					-c <col> # multiple columns separated by commas in same order as of bed format
					-s <separatedBy>
					-r <header> # header to be used in the bed formatted file
					-f <first_line>
					-h <help> 
				'''
			sys.exit(2)
		elif opt in ("-i", "--file"):
			file = arg
		elif opt in ("-c", "--col"):
			col = arg
		elif opt in ("-s", "--separatedBy"):
			sep = arg
		elif opt in ("-f","--first_line"):
			first_line = True
		elif opt in ("-r","--header"):
			header = arg

	return file, col, sep, first_line, header
	

### hash the first file
def PRINT(file,c,sep,first_line, header):
	print header
	for line in open(file,'r'):
		if len(line) > 0:
			if line[0] != '#':
				if first_line == False:
					line = line.strip()
					token = line.split(sep)
					lis = list(col.split(','))
					key = ''
					for i in lis:
						if token[int(i)-1] == 'F':
							value = '+'
						elif token[int(i)-1] == 'R':
							value = '-'
						else:
							value = token[int(i)-1]
						key += value+'\t'
					print key
				first_line = False
	return hash


if __name__ == "__main__":
    
	file, col, sep, first_line, header = options(sys.argv[1:])
	
	### print bed format file
	PRINT(file,col,sep,first_line, header)
