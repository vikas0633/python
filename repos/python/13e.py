### 10b.py for removing additional N region

import sys
st=sys.argv[3]
seq=''
number_of_bases=int(sys.argv[2])
insert=number_of_bases/2
first_line = True
for line in open(sys.argv[1], 'r'):
	line = line.strip()
	try:
		if((line[0] == '>')&(first_line == True)) :
			line=line.strip()
			token=line.split(' ')
			seq_name = token[0][0]+str(st)+'_'+token[0][1:]+'_'+str(insert) ### adding bars to genome name
		if((line[0] == 'A')|(line[0] == 'T')|(line[0] == 'G')|(line[0] == 'C')|(line[0] == 'N')):
			seq += line
		elif(first_line == False):
			
			b_seq = seq[0:(number_of_bases+1)]
			b_token = b_seq.split("NNNN")
			
			c_seq = seq[((number_of_bases) + 1 ): (len(seq)-(number_of_bases) + 1 )]
			a_seq = seq[(len(seq)-(number_of_bases)+1 ):]
			a_token = a_seq.split("NNNN")
			
			f_seq = b_token[len(b_token) - 1] + c_seq + a_token[0]
			print(last_seq_name+','+str(len(b_token[len(b_token) - 1]))+','+str(len(a_token[0]))+','+str(len(c_seq)))
			print(f_seq)
			seq = ''
			seq_name = line
			break
		first_line = False
		last_seq_name = seq_name
	except:
		continue
