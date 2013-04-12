#12ai.py - make a file with unique abundances

import os,sys

### load the genotype name
work_dir=os.getcwd()
lib_path = os.path.abspath(os.getcwd())
script_dir=os.getcwd()+'/'+'scripts'
sys.path.append(lib_path)
from configuration_12 import genotype


def parse_file():
	### load unique read counts
	first_line=True
	abu = {} ### a list to store abundace  
	len_type={} ### store different type of length
	for line in open(infile,'r'):
		line=line.strip()
		token=line.split('\t')
		if first_line==True:
			string = 'Size'
			for name in genotype:
				abu[name]={}
				string += ('\t'+name)
			o.write(string+'\n')
		else:	
			len_type[len(token[0])]=''
			for i in range(len(genotype)):
				if int(token[i+1])!=0:			
					if len(token[0]) in abu[genotype[i]]:
						abu[genotype[i]][len(token[0])] += int(token[i+1])
					else:
						abu[genotype[i]][len(token[0])] = int(token[i+1])		
		first_line=False
	
	for key2 in sorted(len_type.iterkeys()):
		string=''
		for key1 in genotype:
			string += '\t'+str(abu[key1][key2])
		
		o.write(str(key2)+string+'\n') 

	o.close()

if __name__ == "__main__":
	
	
	infile=sys.argv[1]
	o=open(sys.argv[2],'w')
	parse_file()