
# Uses : python 20101206_cluster_analysis.py /Users/vikas0633/Desktop/plant/20101204_gtf/all_compose_sorted.gtf /Users/vikasgupta/Desktop/plant/20101115_sRNA_cluster/all_sRNA_cluster.bed cluster_position_on_genome.xls


import sys, time
start = time.clock()
infile2 = sys.argv[2] ### cluster.bed
count=0
def hash(last_id,id,last_ch,last_st,last_en,ch,st,en):
	if(id!=last_id):
		exon[st,en,ch]='exonic-intergenic overlap'
		inter[last_en,st,ch]='exon-intergenic overlap'
	else:
		exon[st,en,ch]='exon-intron overlap'
		intron[last_en,st,ch]='exon-intron overlap'


	return id,ch,st,en,exon,intron,inter


### hash the compose genome into hash
exon={}
intron={}
inter={}

last_id=''
last_ch=0
last_st=0
last_en=0
### all_compose_sorted.gtf 
for line in open(sys.argv[1],'r'):
	line = line.strip()
	token= line.split('\t')
	if token[2] == 'exon':
		id=token[8]
		ch= token[0]
		st=int(token[3])
		en=int(token[4])
		last_id,last_ch,last_st,last_en,exon,intron,inter=hash(last_id,id,last_ch,last_st,last_en,ch,st,en)
	

	
### cluster.bed
header=True
for line in open(sys.argv[2],'r'):
	line=line.strip()
	flag=0
	if(header==True):
		header=False
		print (line+'\t'+'genomic region')
	else:
		token=line.split('\t')
		ch=token[0]
		st=int(token[1])
		en=int(token[2])
	
		### check where start lies
		for key in exon:
			if((key[0] < en < key[1])&(key[2]==token[0])):
				if(key[0] < st):
					print(line+'\t'+'exonic')
					flag=1
				else:
					print(line+'\t'+exon[key])
					flag=1
			else:
				continue
		for key in intron:
			if((key[0] < en < key[1])&(key[2]==token[0])):
				if(key[0] < st):
					print (line+'\t'+"intronic")
					flag=1
				else:
					print(line+'\t'+intron[key])
					flag=1
			else:
				continue
		for key in inter:
			if((key[0] <= en <= key[1])&(key[2]==token[0])):
				if(key[0] <= st):
					print(line+'\t'+'intergenic')
					flag=1
				else:
					print(line+'\t'+inter[key])
					flag=1
			else:
				continue
		if(flag==0):
			print(line+'\t'+'overlapping')