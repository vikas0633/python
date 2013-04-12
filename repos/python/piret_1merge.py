### merge the columns

import sys
import os.path
if os.path.isfile('*.piret'):
	os.system('rm *.piret')
file_no = 0
line_count = -1
first_line = True
for line in open(sys.argv[1],'r'):
	line = line.strip()
	line_count += 1
	if line_count%111 == 0:
		file_no += 1
		if first_line == False:
			o.close()
		o = open(str(file_no)+'.piret','w')
	if len(line) > 1:	
		token = line.split('\t')
		o.write(token[1]+'\n')
	first_line = False
o.close()
os.system('paste `ls *.piret|sort -n` >merged.txt')
os.system('rm *.piret')
