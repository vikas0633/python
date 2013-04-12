#12ad.py - script for plotting abundaces

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
first_line=True
ab=[] ### store abundace
count=-1
names=[]
for line in open(sys.argv[1],'r'):
	line=line.strip()
	if(first_line==False):
		count += 1
		ab.append([])
		line=line.strip()
		token=line.split('\t')
		for i in range(len(genotype)):
			names.append(genotype[i])
			ab[count].append(int(token[i+1]))
	first_line=False
lib_ab=map(sum, zip(*ab))
max_lib_ab=max(lib_ab)




############################################################
################# Unique reads         #####################
############################################################
first_line=True
ab_uni=[] ### store abundace
count=-1
names=[]
for line in open(sys.argv[2],'r'):
	line=line.strip()
	if(first_line==False):
		count += 1
		ab_uni.append([])
		line=line.strip()
		token=line.split('\t')
		for i in range(len(genotype)):
			names.append(genotype[i])
			ab_uni[count].append(int(token[i+1]))
	first_line=False
lib_ab_uni=map(sum, zip(*ab_uni)) ### unique read counts


### make raw abundace plot
N=len(genotype)
n=num_lines/2
ind = np.arange(N)    # the x locations for the groups
width = 0.35       # the width of the bars: can also be len(x) sequence
#plt.subplots_adjust(wspace=0.4,hspace=0.3)
#plt.subplot(n,2,l)
#print lib_ab,lib_ab_uni
p1=plt.bar(ind,lib_ab,width,color=(38/float(255),163/float(255),91/float(255)))
p2=plt.bar(ind+width,lib_ab_uni,width,color=(169/float(255),237/float(255),44/float(255)))
#plt.xlabel('Libraries',fontsize='small')
matplotlib.rc('ytick', labelsize=18)
plt.title('Filtered Raw Read Counts')
plt.ylabel('Read Counts',fontsize=16)
plt.legend( (p1[0], p2[0]), ('Redundant', 'Non-redundant'),bbox_to_anchor=(1.1, 1.05))
plt.xticks(ind+width, names ,fontsize=7,rotation=30)
plt.savefig(sys.argv[1]+'_'+"Adapter_filtered_total_abundances.png")
plt.clf()

lib_ab=[float(max_lib_ab)/i for i in lib_ab]  ### normalize raw_count

### make plots for librarywise sequence distribution
first_line=True
l=0
abu=[]
abu_uni=[]
size_dic={}

for (line1,line2) in zip(open(sys.argv[1],'r'),open(sys.argv[2],'r')):
	x=[];y=[];names=[]
	if(first_line==False):
		import numpy as np
		import matplotlib.pyplot as plt
		abu.append([])
		abu_uni.append([])
		l += 1
		line1=line1.strip()
		token1=line1.split('\t')
		line2=line2.strip()
		token2=line2.split('\t')
		for i in range(len(genotype)):
			abu[len(abu)-1].append(int(token1[i+1])*lib_ab[i])
			abu_uni[len(abu)-1].append(int(token2[i+1])*lib_ab[i])
			x.append(int(token1[i+1])*lib_ab[i])
			y.append(int(token2[i+1])*(lib_ab[i]))
			names.append(genotype[i])
		fig = plt.figure()
		ax = fig.add_subplot(111)
		N=len(genotype)
		n=num_lines/2
		fig.set_size_inches(12.5,9)
		ind = np.arange(N)    # the x locations for the groups
		width = 0.35       # the width of the bars: can also be len(x) sequence
		x=[item*100/sum(x) for item in x]
		y=[item*100/sum(y) for item in y]
		rects1 = ax.bar(ind,tuple(x),width,color=(1,float(143)/255,float(91)/255))
		rects2 = ax.bar(ind+width,tuple(y),width,color=(1,float(161)/255,float(235)/255))
		#plt.xlabel('Libraries',fontsize='small')
		size_dic[(token1[0])]=''
		matplotlib.rc('ytick', labelsize=18)
		ax.set_ylabel('Reads fraction in % from read length '+ str(token1[0]),fontsize=16)
		plt.xticks(ind+width/4,rotation='30',fontsize=11)
		ax.set_xticklabels((names))
		plt.title('Normalized Read Counts',fontsize=16)
		def autolabel(rects,anno,wd):
			# attach some text labels
			for rect in rects:
				height = rect.get_height()
				plt.text(rect.get_x()+rect.get_width()/2.+wd, 0.5*height, str(anno),
						ha='center', va='bottom')
		
		#autolabel(rects1,'R',0)
		#autolabel(rects2,'U',0)
		
		ax.legend( (rects1[0], rects2[0]), ('Redundant', 'Non-redundant'),bbox_to_anchor=(1.1, 1.05))
		plt.savefig(sys.argv[1]+'_'+str(token1[0])+"_total_abundances.png")
		plt.show()
		plt.clf()
	first_line=False
#plt.title('Unique Size distrubution across libraries')

### change size from dict to list
size=[key for key in sorted(size_dic.iterkeys())]
#print size


##### make librarywise size distribution plots
l=0
for g in range(len(genotype)):
	set=[];	set_uni=[];
	
	for s in range(len(size)):
		set.append(abu[s][g])
		set_uni.append(abu_uni[s][g])
	N=len(size)
	ind = np.arange(N)    # the x locations for the groups
	width = 0.4      # the width of the bars: can also be len(x) sequence
	
	set=[item*100/sum(set) for item in set]
	set_uni=[item*100/sum(set_uni) for item in set_uni]
	p1 = plt.bar(ind+0.5,set,width, color=(67/float(255),146/float(255),1))
	p2 = plt.bar(ind+0.5+width,set_uni,width, color=(109/float(255),212/float(255),1))
	def autolabel(rects,list):
		# attach some text labels
		for rect in rects:
			height = rect.get_height()
			plt.text(rect.get_x()+rect.get_width()/2., 0.5*height, str(round(100*height/float(sum(list)),1))+'%',
					ha='center', va='bottom')
	
	autolabel(p1,set)
	autolabel(p2,set_uni)
	plt.legend( (p1[0], p2[0]), ('Redundant', 'Non-redundant'),bbox_to_anchor=(1.1, 1.05))
	plt.xlabel('Read Length',fontsize=18)
	plt.ylabel('Reads fraction in %',fontsize=18)
	plt.title('Read Length distrubution for '+str(genotype[g]),fontsize=18)
	plt.xticks(ind+width+0.5, size,fontsize=18)
	matplotlib.rc('ytick', labelsize=18)
	plt.savefig(sys.argv[1]+'_'+genotype[g]+'_length_distribution.png')
	plt.clf()

##### make overall size distribution plots

set=[sum(list) for list in abu]
set_uni=[sum(list) for list in abu_uni]
set=[item*100/sum(set) for item in set]
set_uni=[item*100/sum(set_uni) for item in set_uni]
p1 = plt.bar(ind+0.5,set,width, color=(67/float(255),146/float(255),1))
p2 = plt.bar(ind+0.5+width,set_uni,width, color=(109/float(255),212/float(255),1))

autolabel(p1,set)
autolabel(p2,set_uni)
plt.legend( (p1[0], p2[0]), ('Redundant', 'Non-redundant'),bbox_to_anchor=(1.1, 1.05))
plt.xlabel('Read Length',fontsize=18)
plt.ylabel('Reads fraction in %',fontsize=18)
plt.title('Read Length distrubution',fontsize=18)
matplotlib.rc('ytick', labelsize=18)
plt.xticks(ind+width+0.5, size, fontsize=18 )
plt.title('All_data_length_distribution',fontsize=18)
plt.savefig(sys.argv[1]+'_'+'All_data_length_distribution.png')
plt.clf()

