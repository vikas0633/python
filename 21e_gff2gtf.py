#21e_gff2gtf.py - this script converts gff, gff3 format to gtf format - /Users/vikas0633/Desktop/plant/2012_week30

import sys

inFile = open(sys.argv[1],'r')

for line in inFile:
  #skip comment lines that start with the '#' character
  if (line[0] != '#') & (len(line)>1):
    #split line into columns by tab
    data = line.strip().split('\t')
    ID = ''

    #if the feature is a gene 
    if data[2] == "mRNA":
      #get the id
      t_ID = data[-1].split('ID=')[-1].split(';')[0]
      g_ID = data[-1].split('Parent=')[-1].split(';')[0]

        
    ### write mRNA as transcript
    if data[2] == 'mRNA':
      data[2] = 'transcript'
      

    #print out this new GTF line
    if data[2] != 'gene':
      #modify the last column
      data[-1] = 'gene_id "' + g_ID + '"; transcript_id "' + t_ID +'";'

      print '\t'.join(data)
