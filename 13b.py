#### make header of reference sequence suitable for program


import sys

seq=''
first_line=True
for line in open(sys.argv[1],'r'):
	line=line.strip()
	if(line[0]=='>'):
		seq_name=line[1:]
		### check length of seq name
		token=seq_name.split('_')
		if(len(token)<3):
			for i in range(len(token)+1,4):
				seq_name += ('_'+str('test'))
		else:
			seq_name=token[0]
			for j in range(1,3):
				seq_name += ('_'+str(token[j]))
		print(">"+seq_name)
	else:
		seq=line.upper()
		print seq
	first_line=False
			