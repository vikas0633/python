#12s.py - add anotation to the sequences
#nice -n 19 python 12s.py $infile $input.sam 

import sys
annotation=((sys.argv[3]).split('/'))[-1]
### open file with mappings
seq={}
for line in open(sys.argv[2],'r'):
	line = line.strip()
	token = line.split('\t')
	tokens=token[0].split('_')
	key=tokens[1]
	if key in seq:
		seq[key] = (seq[key]+','+token[2]+'_'+token[3])
	else:
		seq[key] = token[2]+'_'+token[3]
		
		
for key in seq:
	token=seq[key].split(',')
	seq[key] = (str(len(token))+';'+str(seq[key]))
		

### open file with sequnces

for line in open(sys.argv[1],'r'):
	line = line.strip()
	token = line.split('\t')
	key=token[0]
	if(line[0]=='#'):
		print ('test_key'+'\t'+'mappings_'+annotation)
	else:
		if key in seq:
			print (key+'\t'+str(seq[key]))
		else:
			print (key+'\t'+'0')