#20_trim_reads.py - script for trimming the fastq reads and quality scores

## Usages: nice -n 19 python 20_trim_reads.py 5'_trim 3'_trim

import sys,os,glob

work_dir = os.getcwd()
path = work_dir+'/'+'output'

### bases to trim from 5'
up =int(sys.argv[1])
### bases to trim from 3'
down =int(sys.argv[2])

### go to output folder and read all the files one by one
for infile in glob.glob(os.path.join(path, '*.fastq')):
    o = open(infile+'_trimmed','w')
    line_no = 0
    for line in open(infile,'r'):
        line_no += 1
        line = line.strip()
        if (line_no%4 == 1) | (line_no%4 == 3):
            o.write(line+'\n')
        else:
            o.write(line[up:len(line)-down]+'\n')
    o.close()