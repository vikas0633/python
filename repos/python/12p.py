#12p.py - /plant/2011_week42/20111019 - for plotting expression data between two sets

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import math
import os,sys
import glob
import re
import time
import datetime
import commands
import linecache
import scipy

work_dir=os.getcwd()
lib_path = os.path.abspath(os.getcwd())
script_dir=os.getcwd()+'/'+'scripts'
sys.path.append(lib_path)
from configuration_12 import genotype
genotype=list(genotype)

def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

### pearson coeff definition
def pearson_def(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)

#output files
o1=open(sys.argv[4],'w') ### values with log odd score
o2=open(sys.argv[5],'w') ### values with relative fraction
#make hash for reference libraries
ref={}
i=0
for line in open(sys.argv[3],'r'):
	line=line.strip()
	if(len(line)>0):
		if(line[0]!='#'):
			i+=1
			token=line.split('\t')
			ref[i]=int(token[1])


#find min in column
x=[]
header=False
for line in open(sys.argv[1],'r'):
	line=line.strip()
	if(len(line)>0):
		if((line[0]=='#')&(header==False)):
			token=line.split('\t')
			l=(len(token)-4)/2
			for i in range(0,l+2):
				x.append([])
				header=True
		elif((line[0]=='#')&(header==True)):
			continue
		else:
			token=line.split('\t')
			start=l+2
			for i in range(start,len(token)):
				if(float(token[i]) > 0):
					x[i-start].append(float(token[i]))			
minimum=[]
for i in range(0,l+2):
	minimum.append(min(x[i]))
		
#make hash for profile libraries
x=[]
y=[]
for line in open(sys.argv[1],'r'):
	line=line.strip()
	out1='';out2='';out='';
	if(len(line)>0):
		if(line[0]=='#'):
			token=line.split('\t')
			l=(len(token)-4)/2
			out += token[0]
			for j in range(l+2,0,-1):
				out += '\t'+str(token[-j])
			out += '\n'
			o1.write(out)
			o2.write(out)
			for i in range(0,l+2):
				x.append([])
				y.append([])
		else:
			token=line.split('\t')
			out1 += token[0]
			out2 += token[0]
			start=l+2
			for i in range(start,len(token)):
				if(float(token[i])==0):
					val=minimum[i-start]
				else:
					val=float(token[i])
				if(int(ref[i+1-start])==0):
					value=1
					norm=val
				else:
					value=float(token[int(ref[i+1-start])+start-1])
					if(value==0):
						value=minimum[i-start]
					norm=((val-value)/(val+value))
				#print ((i+1-start),value,token[i],int(ref[i+1-start]))
				
				x[i-start].append(math.log10(val/value))
				y[i-start].append(norm)
				out1 += '\t'+str(math.log10(val/value))
				out2 += '\t'+str(norm)
			out1 += '\n'
			out2 += '\n'
			o1.write(out1)
			o2.write(out2)
						
			
#print "number of sequences in each library"
#for i in range(0,l+2):
#	print (x[i])
#open file for selective plotting
genotype.append('sum')
genotype.append('score')
line_number=0
o1=open(sys.argv[6],'w')
o2=open(sys.argv[7],'w')
for line in open(sys.argv[2],'r'):
	line=line.strip()
	if(len(line)>0):
		if(line[0]=='#'):
			continue
		else:
			line_number += 1
			token=line.split('\t')
			for i in range(1,len(token)):
				flag=int(token[i])
				if(len(genotype)!=0):
					name=genotype
				else:
					name=range(1,len(token))
				if(flag==1):
					from scipy.stats.stats import *
					plt.plot(x[i-1],x[line_number -1],'.',color='r',markersize=2.5)
					plt.ylabel(str(name[line_number-1])+' normalized expression')
					plt.xlabel(str(name[i-1])+' normalized expression')
					plt.title('Expression distrubution log-score Spearman Corr. Coeff.: '+str(spearmanr(x[line_number -1],x[i-1])[0]),fontsize='small')
					plt.savefig(sys.argv[1]+'_log_'+str(line_number)+'_'+str(i)+".png")
					o1.write(str(name[i-1])+'\tvs\t'+str(name[line_number-1])+'\t'+str(spearmanr(x[line_number -1],x[i-1])[0])+'\n')
					plt.clf()
					plt.plot(y[i-1],y[line_number -1],'.',color='r',markersize=2.5)
					plt.ylabel(str(name[line_number-1])+' normalized expression')
					plt.xlabel(str(name[i-1])+' normalized expression')
					plt.title('Expression distrubution relative fraction Spearman Corr. Coeff.: '+str(spearmanr(y[line_number -1],y[i-1])[0]),fontsize='small')
					plt.savefig(sys.argv[1]+'_'+str(line_number)+'_'+str(i)+".png")
					o2.write(str(name[i-1])+'\tvs\t'+str(name[line_number-1])+'\t'+str(spearmanr(y[line_number -1],y[i-1])[0])+'\n')
					plt.clf()
o1.close()
o2.close()
