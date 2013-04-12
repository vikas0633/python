### check if you have all elements
### Usages: python 10r.py mapping_count_file gap_position_file >output_file

import sys
i=0
line={}
### open file with the mapping counts
for lines in open(sys.argv[1],'r'):
	lines=lines.strip()
	token=lines.split('\t')
	if(len(token)>1):
		i += 1
		if(i==1):
			line[i,token[0]]=lines
		if(i==2):
			line[i,token[0]]=lines
		if(i==3):
			line[i,token[0]]=lines
		if(i==4):
			line[i,token[0]]=lines
			i=0

### open file containing all the gap regions
for lines in open(sys.argv[2],'r'):
	lines=lines.strip()
	key=1,lines
	if key in line:
		for i in range(1,5):
			print(line[i,lines])
	else:
		for i in range(1,5):
			print(lines+'\t'+'None'+'\t'+'0'+'\t'+'0'+'\t'+'500	0	0	0	250	0	0	0	2000	0	0	0	5000	0	0	0')
