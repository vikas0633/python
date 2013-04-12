### 23a_mysql_header.py -/Users/vikas0633/Desktop/script/python/ - script for making headers for mysql tables
### python 23a_mysql_header.py /Users/vikas0633/Downloads/2012_09_18_Mloti_sRNAs_NodEInf_psRNATargetOutput.txt

import os, sys


def is_number(s):
    try:
        float(s)
        return 'FLOAT'
    except ValueError:
        return 'VARCHAR(100)'


infile = sys.argv[1]
print infile

first_line = True
for line in open(infile,'r'):
	token = line.strip().split('\t')
	if first_line == True:
		header = token
	else:
		print 'CREATE TABLE `'+infile.split('/')[-1]+ '` ( '
		
		for i in range(len(token)-1):
			print '`'+str(header[i].replace('"',''))+'` '+is_number(token[i])+' ,'
		i = len(token)-1
		print '`'+str(header[i].replace('"',''))+'` '+is_number(token[i])
		print ' );' 
		break
	first_line = False