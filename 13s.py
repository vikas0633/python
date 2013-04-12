#13s.py -    add distances from 5 prime and 3 prime ends
#nice -n 19 python 13s.py "$input"summary_"$dat"_len $pair_read_2.unmapped.fq.unique_all_rep_aln.sam $pair_read_3.unmapped.fq.unique_all_rep_aln.sam

import sys
### open contig file to calculate contig length
seq=''
key=''
length={}
for line in open(sys.argv[4],'r'):
	line=line.strip()
	if(len(line)>0):
		if(line[0]=='>'):
			if(len(key)!=0):
				length[key]=len(seq)
			key=line[1:]
			seq=''
		else:
			seq += line
length[key]=len(seq)

length['None']=0
#print "contig lenght done"
### open sam file to load data into hash
sam={} ### keys=pos,contigs,orientation,end

for line in open(sys.argv[2],'r'):
	if(len(line)>0):
		line = line.strip()
		token = line.split('\t')
		contig = token[2]
		map = int(token[3])
		tokens = token[0].split('_')
		pos = tokens[0]
		
		### consider contig orientation as forward 
		### distance from 5'
		key=pos,contig,"same","up"
		if key in sam:
			sam[key] += (','+str(map))
		else:
			sam[key] = str(map)
				### distance from 5'
		key=pos,contig,"reverse","up"
		if key in sam:
			sam[key] += (','+str(length[contig]-map))
		else:
			sam[key] = str(length[contig]-map)
		
#print "5 sam done"
for line in open(sys.argv[3],'r'):	
	line = line.strip()
	if(len(line)>0):
		### distance from 3'
		token = line.split('\t')
		contig = token[2]
		map = int(token[3])
		tokens = token[0].split('_')
		pos = tokens[0]
		key=pos,contig,"same","down"
		if key in sam:
			sam[key] += (','+str(length[contig]-map))
		else:
			sam[key] = str(length[contig]-map)
		### consider contig orientation as reverse 
		### distance from 3'
		key=pos,contig,"reverse","down"
		if key in sam:
			sam[key] += (','+str(map))
		else:
			sam[key] = str(map)
# print "3 sam done"
### open final output file and insert the columns
for line in open(sys.argv[1],'r'):
	if(len(line)>0):
		temp=''
		line=line.strip()
		token=line.split('\t')
		contig=token[4]
		if(line[0]=='#'):
			temp = token[0]+'\t'+token[1]+'\t'+token[2]+'\t'+token[3]+'\t'+token[4]+'\t'+token[5]+'\t'
			temp = temp + '5-distance' + '\t'+ '3-distance'
			for i in range(6,len(token)):
				temp = temp + ('\t'+token[i])
			print temp
		else:
			temp = token[0]+'\t'+token[1]+'\t'+token[2]+'\t'+token[3]+'\t'+token[4]+'\t'+token[5]+'\t'
			pos=(10000000000*int(token[2]))+int(token[0])
			try:
				tokens=sam[str(pos),contig,token[7],"up"].split(',')
				avg_up=0
				for i in range(0,len(tokens)):
					avg_up += int(tokens[i])
				avg_up = avg_up/len(tokens)

			except:
				avg_up=0
			try:
				tokens=sam[str(pos),contig,token[7],"down"].split(',')
				avg_down=0
				for i in range(0,len(tokens)):
					avg_down += int(tokens[i])
				avg_down = avg_down/len(tokens)

			except:
				avg_down=0
			temp = temp + str(avg_up) + '\t'+ str(avg_down)
			for i in range(6,len(token)):
				temp = temp + ('\t'+token[i])
			print temp