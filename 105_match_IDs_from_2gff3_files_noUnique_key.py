### 105_match_IDs_from_2gff3_files.py - /Users/vikas0633/Desktop/script/python/ - script will take two gff3 files and print out the corresponding mRNA IDs

import os,sys,getopt, re

### Usage: python ~/script/python/105_match_IDs_from_2gff3_files.py -i $gff3_1 -j $gff3_2 -c 1,3,4,5 -d 1,3,4,5

### main argument to 

def options(argv):
	file1 = ''
	file2 = ''
	col1 = ''
	col2 = ''
	sep='\t'
	first_line = False
	try:
		opts, args = getopt.getopt(argv,"hi:j:c:d:s:",["file1=","file2=","col1=","col2=",'separatedBy='])
	except getopt.GetoptError:
		print '''
			python 100_intersect_columns.py
				-i <file1>
				-j <file2>
				-c <col1> # multiple columns separated by commas
				-d <col2> # multiple columns separated by commas
				-s <separatedBy>
				-f <first_line>
				-h <help> 
			'''
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print '''
					python 100_intersect_columns.py
						-i <file1>
						-j <file2>
						-c <col1> # multiple columns separated by commas
						-d <col2> # multiple columns separated by commas
						-f <first_line>	# Use this option to skip fist_line/header					
						-s <separatedBy>
						-h <help> 
				'''

			sys.exit(2)
		elif opt in ("-i", "--file1"):
			file1 = arg
		elif opt in ("-j", "--file2"):
			file2 = arg
		elif opt in ("-c", "--col1"):
			col1 = arg
		elif opt in ("-d", "--col2"):
			col2 = arg
		elif opt in ("-s", "--separatedBy"):
			sep = arg
		elif opt in ("-f","--first_line"):
			first_line = True
			
	return file1, file2, col1, col2, sep, first_line
	

### hash the first file
def HASH(file1,c1,sep,first_line):
	hash = {}
	count = 0
	lis = []
	for line in open(file1,'r'):
		
		if len(line) > 1:
			if line[0] != '#':
				if first_line == False:
					line = line.strip()
					token = line.split(sep)
					lis = list(col1.split(','))
					key1 = ''
					for i in lis:
						key1 += '-'+token[int(i)-1]
					if (token[2]=='mRNA') and (token[1]!='rRNA') and (token[1]!='tRNA'):
						match = re.search(r'ID=.+;',line)
						if match:
							match1 = match.group().split(';')[0].replace(';','').replace('ID=','')
							key1 += '-'+match1.split('.')[-1]
							for line in open(file2,'r'):
								if len(line) > 1:
									if line[0] != '#':
										if first_line == False:
											line = line.strip()
											token = line.split(sep)
											lis = list(col2.split(','))
											key2 = ''
											for i in lis:
												key2 += '-'+token[int(i)-1]
											if (token[2]=='mRNA') and (token[1]!='rRNA') and (token[1]!='tRNA'):
												match = re.search(r'ID=.+;',line)
												if match:
													match2 = match.group().split(';')[0].replace(';','').replace('ID=','')
													key2 += '-'+match2.split('.')[-1] 
												if match2 not in lis: 
													if key1==key2:
														print match1+'\t'+match2
														lis.append(match2)
														break
										first_line = False



if __name__ == "__main__":
    
	file1, file2, col1, col2, sep, first_line = options(sys.argv[1:])
	
	### hash the first file
	hash = HASH(file1,col1,sep,first_line)
	
	### parse the second file
	PARSE(file2,col2,sep,first_line, hash)
