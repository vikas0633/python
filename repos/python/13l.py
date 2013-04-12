####### summary script for rep elements

import sys
sum1={}
sum2={}
last_ori=1;ori=0
last_gap_pos=0;gap_pos=0
first_line=True
score={}
lb_size={1:500,2:250,3:2000,4:5000}
reward={1:100,2:200,3:1,4:1}### rewarding for smaller size library
for line in open(sys.argv[1],'r'):
	line=line.strip()
	already_seen=False
	count=0
	if((len(line)>0) & (line[0:3]!='END')):
		token=line.split('\t')
		tokens=token[0].split('_')
		gap_pos=int(tokens[0])
		rep_name=token[1]
		rep_count=float(token[2])
		if(float(rep_count)==0):
			rep_count=0
		ori=float(token[3])
		#print(last_gap_pos,gap_pos)
		if((first_line==False)&((last_ori!=ori)|(last_gap_pos!=gap_pos))): ### only for ori =1, stands for both direction reads
			key=max(sum1,key=sum1.get)
			score[int(key[0])]=str(key[1]),sum1[key]
			sum1={}
		sum1[gap_pos,rep_name,rep_count]=rep_count ###accounts for number of hanging reads
		for i in range(1,5):
			if((float(token[4*i+2])==0)):
				token[4*i+2]=50000
			if((float(token[4*i+3])==0)):
					token[4*i+3]=1
			if(already_seen==True):
				count += 100*float(token[4*i+1])*(1+(1/(float(token[4*i+2])+float(token[4*i+3]))))*reward[i]
			if(already_seen==False):
				count += float(token[4*i+1])*(1+(1/(float(token[4*i+2])+float(token[4*i+3]))))*reward[i]
			if((float(token[4*i+1])!=0)):
				already_seen=True
			#count += float(token[4*i+1])*float(token[4*i+1])*(1+(1/(float(token[4*i+2])+float(token[4*i+3]))))		
			#count += float(token[4*i+1])*(float(token[4*i+2])*float(token[4*i+3]))/(lb_size[i]*lb_size[i])
			#count += float(token[4*i+1])
			#count += float(token[4*i+1])/(float(token[4*i+2])+float(token[4*i+3]))
		sum1[gap_pos,rep_name,rep_count]=float(count)
		last_gap_pos=gap_pos
		last_ori=ori
		first_line=False
	if((line[0:3]=='END')): ###new change
		key=max(sum1,key=sum1.get)
		score[int(key[0])]=str(key[1]),sum1[key]
		
last_ori=1;
last_gap_pos=0;gap_pos=0
first_line=True

### use second file and take out better orientation	
for line in open(sys.argv[2],'r'):
	line=line.strip()
	count=0
	already_seen=False
	if((len(line)>0) & (line[0:3]!='END')):
		token=line.split('\t')
		tokens=token[0].split('_')
		gap_pos=int(tokens[0])
		rep_name=token[1]
		rep_count=float(token[2])
		if(float(rep_count)==0):
			rep_count=0
		ori=float(token[3])
		if((first_line==False)&((last_ori!=ori)|(last_gap_pos!=gap_pos))): ### only for ori =1, stands for both direction reads
			key=max(sum2,key=sum2.get)
			#print(float(sum2[key]),float(score[key[0]][1]))
			if(float(sum2[key])>=float(score[key[0]][1])):
				print (str(key[0])+'\t'+str(key[1])+'\t'+str(float(sum2[key]))+'\t'+'reverse')
			else:
				print (str(key[0])+'\t'+score[key[0]][0]+'\t'+str(float(score[key[0]][1]))+'\t'+'same')
			sum2={}
		sum2[gap_pos,rep_name,rep_count]=rep_count ###accounts for number of hanging reads
		for i in range(1,5):
			if((float(token[4*i+2])==0)):
				token[4*i+2]=50000
			if((float(token[4*i+3])==0)):
					token[4*i+3]=1
			"""
			sum2[gap_pos,rep_name,rep_count] += (float(token[3*i+1])/(float(token[3*i+2])*float(token[3*i+3])))
		sum2[gap_pos,rep_name,rep_count] += (float(rep_count)/(50000*50000))
		#print(gap_pos,rep_name,sum2[gap_pos,rep_name,rep_count])
		#sum2[gap_pos,rep_name,rep_count] += ((20000*float(token[3*i+1]))-float(token[3*i+2])*float(token[3*i+1])--float(token[3*i+3])*float(token[3*i+1]))
			"""
			if(already_seen==True):
				count += 100*float(token[4*i+1])*(1+(1/(float(token[4*i+2])+float(token[4*i+3]))))*reward[i]
			if(already_seen==False):
				count += float(token[4*i+1])*(1+(1/(float(token[4*i+2])+float(token[4*i+3]))))*reward[i]
			if((float(token[4*i+1])!=0)):
				already_seen=True
			#count += float(token[4*i+1])*(float(token[4*i+2])*float(token[4*i+3]))/(lb_size[i]*lb_size[i])
			#count += float(token[4*i+1])
			#count += float(token[4*i+1])/(float(token[4*i+2])+float(token[4*i+3]))
		sum2[gap_pos,rep_name,rep_count]=float(count)
		last_gap_pos=gap_pos
		last_ori=ori
		first_line=False	
	if((line[0:3]=='END')): ###new change
		key=max(sum2,key=sum2.get)
		#print(float(sum2[key]),float(score[key[0]][1]))
		if(float(sum2[key])>=float(score[key[0]][1])):
			print (str(key[0])+'\t'+str(key[1])+'\t'+str(float(sum2[key]))+'\t'+'reverse')
		else:
			print (str(key[0])+'\t'+score[key[0]][0]+'\t'+str(float(score[key[0]][1]))+'\t'+'same')

