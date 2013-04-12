import os,sys
import math

def is_float(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False

def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

### pearson coeff definition
#def pearson_def(x, y):
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
    
### pearson rank coeff definition 
from scipy.stats.stats import *


if __name__ == "__main__":
	store=[]
	first_time=True
	for line in open(sys.argv[1],'r'):
		line =line.strip()
		## header
		if ((line[0]=='#')&(first_time==True)):
			header=line.split('\t')
			for i in range(0,len(header)):
				store.append([])
			first_time=False
		#""" works on normal way
		#else:
		#	token=line.split('\t')
		#	for i in range(0,len(token)):
		#		if is_float(token[i]):
		#			store[i].append(float(token[i]))
		#"""			
		
		### for columns when we need to take reference into account
		else:
			token=line.split('\t')
			for i in range(0,len(token)):
				if(i>1):
					if is_float(token[i]):
						if(float(token[i-1])==0):
							token[i-1]=1
						store[i].append(float(token[i])/float(token[i-1]))
		
	print "spearman correlation co-efficeint is:"
	for i in range(0,len(store)):
		for j in range(0,len(store)): 
			if((len(store[i]) >0) & (len(store[j]) >0)):
				print header[i] + " vs " + header[j] + " : "+str(spearmanr(store[i],store[j])[0])