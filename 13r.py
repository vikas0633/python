#13r.py -  	add length of the contigs
import sys
### open contig file to calculate contig length
seq=''
key=''
length={}
for line in open(sys.argv[1],'r'):
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

#open file where length column to be added
for line in open(sys.argv[2],'r'):
	line=line.strip()
	temp_line=''
	if(len(line)>0):
		if(line[0]=='#'):
			token=line.split('\t')
			temp_line=str(token[0])+'\t'+str(token[1])+'\t'+str(token[2])+'\t'+str(token[3])+'\t'+str(token[4])+'\t'+"length_of_contig"+'\t'+str(token[5])+'\t'+str(token[6])+'\t'+str(token[7])+'\t'+"length_of_contig"+'\t'+str(token[8])+'\t'+str(token[9])+'\t'+str(token[10])+'\t'+"length_of_contig"+'\t'+str(token[11])+'\t'+str(token[12])
			print temp_line
		else:
			token=line.split('\t')
			temp_line=str(token[0])+'\t'+str(token[1])+'\t'+str(token[2])+'\t'+str(token[3])+'\t'+str(token[4])+'\t'+str(length[token[4]])+'\t'+str(token[5])+'\t'+str(token[6])+'\t'+str(token[7])+'\t'+str(length[token[7]])+'\t'+str(token[8])+'\t'+str(token[9])+'\t'+str(token[10])+'\t'+str(length[token[10]])+'\t'+str(token[11])+'\t'+str(token[12])
			print temp_line