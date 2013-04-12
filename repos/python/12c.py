### date: 20110506
### Auther: Vikas Gupta
### Make_fasta_v1.0

### Input: a file which contains sequence each line
### Output: fasta file containing sequences and line number as header
### Uses: 20110506_make_fasta.py Input Output

import sys

infile = sys.argv[1]
outfile = sys.argv[2]

### open input file
f = open(infile,'r')
### open output file
o = open(outfile,'w')

### counter for line number
line_no = 0
first_line=True
### go through each line
for line in f.readlines():
	line_no += 1
	### ignore header line
	if(first_line==False):
		line = line.strip()
		token = line.split()
		### checking through each token for sequence
		if ((token[0][0]==('A'))|(token[0][0]==('T'))|(token[0][0]==('G'))|(token[0][0]==('C'))|(token[0][0]==('U'))):
			header = '>'+ str(line_no) +'_' +str(token[0]) + '_' + str(len(token[0]))
			o.write(header)
			o.write('\n')
			o.write(token[0])
			o.write('\n')
	first_line=False
f.close() ### closing infile
o.close() ### closing outfile
