### same as 20110914_10c*.py making unmapped read files
import sys



### make hash from sam file rather than fq file may be a bit slower but memory efficient
o = open(sys.argv[3],'w') ### output sam file
line_sam={}
line_genome_id={}
line_map={}
line_upstream={}
for lines in open(sys.argv[1],'r'):
	lines = lines.strip()
	tokens = lines.split('\t')
	token = tokens[0].split('/')
	if(int(token[1])==1):
		seq_id = token[0][0:]+'/'+'2'
	if(int(token[1])==2):
		seq_id = token[0][0:]+'/'+'1'
	token1 = tokens[2].split(',')
	up_stream=int(token1[1])
	o.write(lines);o.write('\n')
	genome_id=tokens[2]+','+str(tokens[3])+','
	line_sam[seq_id]=1
	line_genome_id[seq_id]=genome_id
	line_map[seq_id]=int(tokens[3])
	line_upstream[seq_id]=up_stream
		
o2= open(sys.argv[6],'w') ### unique unmapped sequences
nb= int(sys.argv[4])
ori=int(sys.argv[7]) ###ori
i=0
first_line=True
last_key=''
for line in open(sys.argv[2],'r'):
	line=line.strip()
	tokens = line.split('\t')
	#token = tokens[0].split('/')
	i += 1
	if(line[0]=='@'):
		line1={}
		line2={}
		line3={}
		line4={}
		i=1
		key=tokens[0][1:]
		line1[key]=line ### store first line
	if(line[0]!='@'):
		if(i==2):				### store second line
			line2[key]=line
		if(i==3):				### store third line
			line3[key]=line
		if(i==4):				### store fourth line
			line4[key]=line
	if((first_line==False)&(line[0]=='@')):
		if(ori==1):
			if last_key in line_sam:
				id=line_genome_id[last_key].replace(':','_')
				o2.write(last_line1[last_key][0]+id+'x'+last_line1[last_key][1:]);o2.write('\n');
				print(last_line1[last_key][0]+id+'x'+last_line1[last_key][1:])
				o2.write(last_line2[last_key]);o2.write('\n');
				print(last_line2[last_key])
				o2.write(last_line3[last_key]);o2.write('\n');
				print(last_line3[last_key])
				o2.write(last_line4[last_key]);o2.write('\n');
				print(last_line4[last_key])
		if(ori==2):
			if last_key in line_sam:
				if(line_map[last_key] <= line_upstream[last_key]):
					id=line_genome_id[last_key].replace(':','_')
					o2.write(last_line1[last_key][0]+id+'x'+last_line1[last_key][1:]);o2.write('\n');
					print(last_line1[last_key][0]+id+'x'+last_line1[last_key][1:])
					o2.write(last_line2[last_key]);o2.write('\n');
					print(last_line2[last_key])
					o2.write(last_line3[last_key]);o2.write('\n');
					print(last_line3[last_key])
					o2.write(last_line4[last_key]);o2.write('\n');
					print(last_line4[last_key])
		if(ori==3):
			if last_key in line_sam:
				if(line_map[last_key] > line_upstream[last_key]):
					id=line_genome_id[last_key].replace(':','_')
					o2.write(last_line1[last_key][0]+id+'x'+last_line1[last_key][1:]);o2.write('\n');
					print(last_line1[last_key][0]+id+'x'+last_line1[last_key][1:])
					o2.write(last_line2[last_key]);o2.write('\n');
					print(last_line2[last_key])
					o2.write(last_line3[last_key]);o2.write('\n');
					print(last_line3[last_key])
					o2.write(last_line4[last_key]);o2.write('\n');
					print(last_line4[last_key])
	if(line[0]=='@'):
		last_key=key
	last_line1=line1
	last_line2=line2
	last_line3=line3
	last_line4=line4
	first_line=False
	#print(last_line1,last_line2,last_line3,last_line4)
o.close()