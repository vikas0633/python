#9c.py 2011_gap_filling/2011_week33/ for taking out sequence where only one read is mapped

import sys

mate_pair=sys.argv[5]
number_of_bases = int(sys.argv[4])
o2=open(sys.argv[6],'w')

last_token = {};last_token[0] = '';last_token[4] = '';last_lines=''
read_file = sys.argv[2]

### look for the reads and give back the fastq lines
def make_fastq(seq_id_hash,read_file):
	flag = False
	for lines in open(read_file,'r'):
		lines = lines.strip()
		tokens = lines.split('\t')
		token = tokens[0].split('/')
		key = token[0][5:]
		if key in seq_id_hash:
			print (lines)
			flag = True
		else:
			if((flag == True)&(lines[0] != '@')):
				print (lines)
			else:
				flag = False
		
def make_fastq_unique(seq_id_hash_unique,read_file,o2):
	flag = False
	for lines in open(read_file,'r'):
		lines = lines.strip()
		tokens = lines.split('\t')
		key = tokens[0][5:]
		if key in seq_id_hash_unique:
			o2.write(lines);o2.write('\n');
			
			flag = True
		else:
			if((flag == True)&(lines[0] != '@')):
				o2.write(lines);o2.write('\n');
			else:
				flag = False

	o2.close()					
	### convert to binary number

def Denary2Binary(n):
	'''convert denary integer n to binary string bStr'''
	bStr = ''
	if n < 0:  raise ValueError, "must be a positive integer"
	if n == 0: return '0'
	while n > 0:
		bStr = str(n % 2) + bStr
		n = n >> 1
	return bStr

def int2bin(n, count=24):
	"""returns the binary of integer n, using count number of digits"""
	return "".join([str((n >> y) & 1) for y in range(count-1, -1, -1)])

### fill unmapped reads in hash 

from collections import defaultdict
seq_id_hash = defaultdict(int)
seq_id_hash_unique = defaultdict(int)
o = open(sys.argv[3],'w') ### output sam file
for lines in open(sys.argv[1],'r'):
	lines = lines.strip()
	tokens = lines.split('\t')
	if(tokens[0]==last_token[0]):
		bit_flag = '0000000'+str(Denary2Binary(int(tokens[1])))
		bit_flag_last = '0000000'+str(Denary2Binary(int(last_token[1])))


		#### if mate pair data
		if(mate_pair=='T'):
		
			### if one of the read is unmapped
			if((int(bit_flag[len(bit_flag)-3]) == 1)|(int(bit_flag[len(bit_flag)-4]) == 1)):
			### if first read is mapped
				if(int(bit_flag[len(bit_flag)-4]) == 0):
				### if read is forward
					if(int(bit_flag[len(bit_flag)-6]) == 0):
						if(int(tokens[3])< number_of_bases):
							o.write(last_lines);o.write('\n');o.write(lines);o.write('\n')
							seq_id = tokens[0][4:]
							seq_id_hash[seq_id] = 1
							seq_id_unique = tokens[0][4:]+'/2'
							seq_id_hash_unique[seq_id_unique]=1
					else:
						if(int(tokens[3])> number_of_bases):
							o.write(last_lines);o.write('\n');o.write(lines);o.write('\n')
							seq_id = tokens[0][4:]
							seq_id_hash[seq_id] = 1
							seq_id_unique = tokens[0][4:]+'/2'
							seq_id_hash_unique[seq_id_unique]=1
 

			### if second read is mapped
				else:
					### if read is forward
					if(int(bit_flag_last[len(bit_flag_last)-6]) == 0):
						if(int(tokens[3])> number_of_bases):
							o.write(last_lines);o.write('\n');o.write(lines);o.write('\n')
							seq_id = tokens[0][4:]
							seq_id_hash[seq_id] = 1
							seq_id_unique = tokens[0][4:]+'/1'
							seq_id_hash_unique[seq_id_unique]=1

					else:
						if(int(tokens[3])< number_of_bases):
							o.write(last_lines);o.write('\n');o.write(lines);o.write('\n')
							seq_id = tokens[0][4:]
							seq_id_hash[seq_id] = 1
							seq_id_unique = tokens[0][4:]+'/1'
							seq_id_hash_unique[seq_id_unique]=1
							
		### if end pair data
		else:
			### if one of the read is unmapped
			if((int(bit_flag[len(bit_flag)-3]) == 1)|(int(bit_flag[len(bit_flag)-4]) == 1)):
				### if first read is mapped
				if(int(bit_flag[len(bit_flag)-4]) == 0):
					### if read is forward
					if(int(bit_flag[len(bit_flag)-6]) == 0):
						if(int(tokens[3])< number_of_bases):
							o.write(last_lines);o.write('\n');o.write(lines);o.write('\n')
							seq_id = tokens[0][4:]
							seq_id_hash[seq_id] = 1
							seq_id_unique = tokens[0][4:]+'/2'
							seq_id_hash_unique[seq_id_unique]=1
					else:
						if(int(tokens[3])> number_of_bases):
							o.write(last_lines);o.write('\n');o.write(lines);o.write('\n')
							seq_id = tokens[0][4:]
							seq_id_hash[seq_id] = 1
							seq_id_unique = tokens[0][4:]+'/2'
							seq_id_hash_unique[seq_id_unique]=1
			
			
				### if second read is mapped
				else:
					### if read is forward
					if(int(bit_flag_last[len(bit_flag_last)-6]) == 0):
						if(int(tokens[3])< number_of_bases):
							o.write(last_lines);o.write('\n');o.write(lines);o.write('\n')
							seq_id = tokens[0][4:]
							seq_id_hash[seq_id] = 1
							seq_id_unique = tokens[0][4:]+'/1'
							seq_id_hash_unique[seq_id_unique]=1

					else:
						if(int(tokens[3])> number_of_bases):
							o.write(last_lines);o.write('\n');o.write(lines);o.write('\n')
							seq_id = tokens[0][4:]
							seq_id_hash[seq_id] = 1
							seq_id_unique = tokens[0][4:]+'/1'
							seq_id_hash_unique[seq_id_unique]=1
	last_token = tokens
	last_lines = lines
o.close()
		

make_fastq(seq_id_hash,read_file)
make_fastq_unique(seq_id_hash_unique,read_file,o2)