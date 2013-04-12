#12am.py - script for plotting abundaces before adapter filtering

import sys
import os,sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


work_dir=os.getcwd()
lib_path = os.path.abspath(os.getcwd())
script_dir=os.getcwd()+'/'+'scripts'
sys.path.append(lib_path)
from configuration_12 import genotype

### count number of line 
num_lines = sum(1 for line in open(sys.argv[1]))


############################################################
################# Raw Data          ########################
############################################################
### calculate normalization factor
ab=[] ### store abundace
first_line=True
count=-1
names=[]
for line in open(sys.argv[1],'r'):
	if first_line==False:
		line=line.strip()
		count += 1
		line=line.strip()
		token=line.split('\t')
		names.append(genotype[count])
		ab.append(int(token[1]))
	first_line=False



############################################################
################# Unique reads         #####################
############################################################
ab_uni=[] ### store abundace
count=-1
names=[]
first_line=True
for line in open(sys.argv[1],'r'):
	if first_line==False:
		line=line.strip()
		count += 1
		line=line.strip()
		token=line.split('\t')
		names.append(genotype[count])
		ab_uni.append(int(token[2]))
	first_line=False
### make raw abundace plot
N=len(genotype)
n=num_lines/2
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence
#plt.subplots_adjust(wspace=0.4,hspace=0.3)
#plt.subplot(n,2,l)
#print lib_ab,lib_ab_uni
p1=plt.bar(ind,ab,width,color=(38/float(255),163/float(255),91/float(255)))
p2=plt.bar(ind+width,ab_uni,width,color=(169/float(255),237/float(255),44/float(255)))
#plt.xlabel('Libraries',fontsize='small')
plt.title('Raw Read Counts')
plt.ylabel('Read Counts',fontsize=7)
plt.legend( (p1[0], p2[0]), ('Redundant', 'Non-redundant'),bbox_to_anchor=(1.1, 1.05))
plt.xticks(ind+width, names ,fontsize=7,rotation=30)
plt.savefig("Total_library_abundances_no_adapter_filtering.png")
plt.clf()
