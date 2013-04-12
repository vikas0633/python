### for taking out cluster information
#nice -n 19 fastacmd -d all_lj.tfa -p F -s chr2 -L 29632331,29632562

import sys

last_start = 0
infile = open('p0.001_sRNA21nt_out.txt','r')
for line in infile.readlines():
        line = line.strip()
        token = line.split('\t')
        chromosome = token[0]
        start = int(token[1])
        if( start !=  last_start):
                cluster_file = open(sys.argv[1],'r')
                for lines in cluster_file.readlines():
                        lines = lines.strip()
                        tokens = lines.split('\t')
                        if (chromosome == tokens[0]):
                                if( (int(tokens[1]) <= int(start) < int(tokens[2])) & (500> (int(tokens[2])-int(tokens[1])) >150) ): ### range of Tasi-cluster 150nt to 500nt
                                        out = str('nice -n 19 fastacmd -d '+sys.argv[2]+' -p F -s '+str(tokens[0])+ ' -L '+str(tokens[1])+','+str(tokens[2]))
                                        print (out)
        last_start = start
