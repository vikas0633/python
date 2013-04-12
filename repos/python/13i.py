#10q.py - calculating insert_size/distances
# input line should look like
#119270000000009_lcl|BM2600b_test_test_10927-17077_500,6151,886,xHWIEAS210R_0007:8:1:11225:8918#CAAAAT/1 16      LjRE14  4158
import sys,numpy,math
import operator

last_ref=0;insert={};ref=0
count={};
first_line=True
rep={}
ori=sys.argv[4]
strand=int(sys.argv[5])
### open rep file for calculating length of elements
for line in open(sys.argv[3],'r'):
	line=line.strip()
	if(len(line) >0):
		if((line[0]=='>')&(first_line==False)):
			rep[seq_name]=len(seq)
		if(line[0]=='>'):
			seq=''
			seq_name=line[1:]
		else:
			seq += line
	first_line=False
rep[seq_name]=len(seq)

first_line=True	
o=open(sys.argv[2],'w')
for line in open(sys.argv[1],'r'):
	line=line.strip()
	if(len(line) > 0):
		if(line[0:3]!='END'):
			line=line.strip()
			token=line.split()
			rep_name=token[1]
			rep_pos=int(token[2])
			token1=token[0].split('_')
			token2=token[0].split(',')
			token3=token2[0].split('_')
			ref=int(token1[0])
			ref_length=int(token2[1])
			gap_length=int(token2[3])
			ref_pos=int(token2[4])
			expected_insert=int(token3[-1])
	
		if(((last_ref!=ref)&(first_line==False))|(line[0:3]=='END')):
			s = sorted(count.iteritems(), key=operator.itemgetter(1)) ### s is the key for sorted values
			for i in range(1,5):
				o.write('\n')
				o.write(str(last_ref));o.write('\t')
				if(len(s)>3):
					key=s[-i][0]
					o.write(str(key));o.write('\t')
					o.write(str(count[key]));o.write('\t')
					o.write(str(ori));o.write('\t')
					for l in range(1,5):
						o.write(str(expected_ins[l]));o.write('\t')
						list=[]
						key1=l,key
						if key1 in insert:
							token=str(insert[key1]).split('-')
							for t in range(0,len(token)):
								if(float(token[t])==0):
									continue
								else:
									x=math.fabs((float(token[t]) - float(expected_ins[l]))/float(expected_ins[l]))
									list.append(x)
							if(len(list)==0):
								list=[0]
								o.write( str('0')+'\t'+str((math.fabs(numpy.average(list))))+'\t'+str((numpy.std(list)))+'\t')
								list=[]
							if(len(list)>0):
								o.write(str(len(list))+'\t'+str((math.fabs(numpy.average(list))))+'\t'+str((numpy.std(list)))+'\t')
						else:
							list=[0]
							o.write( str('0')+'\t'+str((math.fabs(numpy.average(list))))+'\t'+str((numpy.std(list)))+'\t')
						
			count={};
			insert={};			
		if(line[0:3]!='END'):
			for j in range(0,5):
				for k in range(0,5):
					count["None"+str(k)]=0
					insert[j,"None"+str(k)]=0
			if(strand==0):
				if(ref_pos > int(token2[1])):
					insert_found=(rep[rep_name])-rep_pos+ ref_pos - ref_length -gap_length
				else:
					insert_found=ref_length - ref_pos + rep_pos
			else:
				if(ref_pos > int(token2[1])):
					insert_found=rep_pos+ ref_pos - ref_length -gap_length
				else:
					insert_found=ref_length - ref_pos + (rep[rep_name])-rep_pos
			key=rep_name
			expected_ins={1:500,2:250,3:2000,4:5000}
			for m in range(1,5):
				if(expected_ins[m]==expected_insert):
					key1=m,key
					if key1 in insert:
						insert[key1] += str("-"+str(insert_found))
					else:
						insert[key1]=str(insert_found)
						
					if key in count:
						count[key] += 1
					else:
						count[key] = 1
		last_ref=ref
		first_line=False
