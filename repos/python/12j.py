### for calculating size distribution

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
fig = plt.figure()
fig.set_size_inches(12.5,9)

seq_len={}
for line in open(sys.argv[1],'r'):
	line=line.strip()
	if(len(line)>0):
		if(line[0]!='#'):
			token=line.split('\t')
			seq=token[0]
			length=len(seq)
			if length in seq_len:
				seq_len[length] += 1
			else:
				seq_len[length] = 1
print ("size\t read_count")
for key in seq_len:
	print str(key)+'\t'+str(seq_len[key])
	

count=0
x=[]
y=[]
for key in seq_len:
	count = seq_len[key]
	x.append(count)
	y.append(key)
N=len(seq_len)
ind = np.arange(N)    # the x locations for the groups
width = 0.5       # the width of the bars: can also be len(x) sequence

p1 = plt.bar(ind+0.5,x,width, color=(32/float(255),142/float(255),250/float(255)))
def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 0.5*height, str(round(100*height/float(sum(x)),2))+'%',
                ha='center', va='bottom')

autolabel(p1)
plt.xlabel('Sequence Length')
plt.ylabel('Sequence Count')
plt.title('Sequence Size distrubution')
plt.xticks(ind+width/2.+0.5, y )
plt.savefig(sys.argv[1]+".png")