### python script for reaplcing libraries header

import sys,os
import numpy as np
##import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt


work_dir=os.getcwd()
lib_path = os.path.abspath(os.getcwd())
script_dir=os.getcwd()+'/'+'scripts'
sys.path.append(lib_path)
from configuration_12 import genotype


### open profile file in write mode
o=open(sys.argv[2],'w')

### open the alias file for profile check.txt
first_line=True
for line in open(sys.argv[1],'r'):
	import re
	if(first_line==True):
		for i in range(len(genotype)):
			#line=re.sub('\\blib-'+str(i+1)+'\\b',genotype[i],line)
			line=re.sub('lib-'+str(i+1)+' ',genotype[i],line)
			line=re.sub('lib-'+str(i+1)+'\t',genotype[i]+'\t',line)
			line=re.sub('lib-'+str(i+1)+'_',genotype[i]+'_',line)
		o.write(line)
	else:
		o.write(line)
	first_line=False
o.close()