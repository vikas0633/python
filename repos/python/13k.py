###10t.py - 2011_gap_filling/2011_week39/20110928 ### make summary output for elements
### Usages: nice -n 19 python 10t.py rep_pos "$dat"_summary.txt.1 "$dat"_summary.txt.2 "$dat"_summary.txt.3

import sys

summary={}
for i in range(2,5):
	for line in open(sys.argv[i],'r'):
		line=line.strip()
		token=line.split('\t')
		pos=int(token[0])
		summary[i,pos]=token[1]+'\t'+str(round(float(token[2])/100,2))+'\t'+token[3]


### print header

print ('#count'+'\t'+'Chromosome'+'\t'+'gap_start'+'\t'+'gap_end'+'\t'+'5-3_best_hit'+'\t'+'Score'+'\t'+'orientation'+'\t'+'5_best_hit'+'\t'+'Score'+'\t'+'orientation'+'\t'+'3_best_hit'+'\t'+'Score'+'\t'+'orientation')
### open rep_pos file
for line in open(sys.argv[1],'r'):
	line=line.strip()
	token=line.split('\t')
	count=int(token[0])
	chro=(token[1]).replace('_test','')
	start=int(token[2])
	en=int(token[3])
	st=10000000000*start
	key=st+count
	string=str(str(count)+'\t'+str(chro)+'\t'+str(start)+'\t'+str(en))
	for i in range(2,5):
		string +=  ('\t'+summary[i,key])
	print (string)	