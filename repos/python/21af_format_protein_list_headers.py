'''
Created on Feb 25, 2013

@author: vgupta
'''

import sys

### 12af_format_protein_list_headers.py - script to get the corresponding headers between corrected and real fasta file
###

### hash the modified headers
file_modified='20121227_conserved_proteins_mRNA_seq.fasta_simple_headers'

header_mod = {}
count = 0 ### hash the headers from number 1 to last

for line in open(file_modified,'r'):
    line = line.strip()
    if len(line) > 1:
        if line[0]=='>':
            count +=1 
            header_mod[count] = line

### hash the original headers
file_original='20121227_conserved_proteins_mRNA_seq.fasta'

header_ori = {}
count = 0 ### hash the headers from number 1 to last

for line in open(file_original,'r'):
    line = line.strip()
    if len(line) > 1:
        if line[0]=='>':
            count +=1 
            header_ori[header_mod[count]] = line.split(' ')[0]
        
### modify the cDNA file
cDNA=sys.argv[1]

for line in open(cDNA,'r'):
    line = line.strip()
    if len(line) > 1:
        if line[0] == '>':
            line = line.replace(' ','')
            token = line.split('.')
            if token[0] in header_ori:
                print (header_ori[token[0]])+'.'+token[1]+'.'+token[2]
                flag = True
            else:
                flag = False
        else:
            if flag == True:
                print line
