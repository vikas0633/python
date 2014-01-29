#-----------------------------------------------------------+
#                                                           |
# 30b_degradome_site_relative2_smallRNAs.py - script for counting small Reads relative to cleavage site      |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# Short script to convert and copy the wheat BACs           |
# Run this in the parent dir that the HEX* dirs exist       |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/100b_fasta2flat.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re, time


### global variables
global ifile, time_start
global ifile, dfile, site_cov, flanking, min_abund, max_abund, base_range, chro_list
start_time = time.time()

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')

### return time eclapsed
def PrinteclapsedTime():
    diff = time.time() - start_time
    minutes, seconds = int(diff)/60, diff % 60
    print('Time taken Min:Sec ==> ' + str(minutes) + ':' + str(round(seconds,2)))


### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 30b_degradome_site_relative2_smallRNAs.py
                                        -i <ifile> ## file with first 4 columns as <#Sequence> <read_size> <sum_norm> <lj_r30.fa>
                                        -a <dfile> ### file with degradome sites
                                        -s <site_cov> ### minimum total abudance on a site flanking region [default: 100]
                                        -f <flanking> ### number of flanking bases around the cleavage sie [default: 25]
                                        -m <min_abund> ### minimum abundance cutoff for a read to be considered [default: 1]
                                        -n <max_abund> ### maximum abundance cutoff for a read to be considered [default: 1000000000 ]
                                        -b <base_range> ### only consider the bases in range from X-Y, i.e., 9-12 [default: 1-len(seq)]
                                        -c <chro_list> ### a list of chromosomes to be considered [default: chr0, chr1,chr2,chr3,chr4,chr5,chr6,chloro,mito]
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile, dfile, site_cov, flanking, min_abund, max_abund, base_range, chro_list
    ifile = ''
    dfile = ''
    site_cov = 100
    flanking = 25
    min_abund = 1
    max_abund = 1000000000
    base_range = -1
    chro_list=['chr0','chr1','chr2','chr3','chr4','chr5','chr6','chloro','mito']
    try:
        opts, args = getopt.getopt(argv,"hi:a:s:f:m:n:b:c:",["ifile=","dfile=","site_cov=","flanking=","min_abund=","max_abund=","base_range=","chro_list="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-a", "--dfile"):
            dfile = arg
        elif opt in ("-s", "--site_cov"):
            site_cov = arg
        elif opt in ("-f", "--flanking"):
            flanking = int(arg)
        elif opt in ("-m", "--min_abund"):
            min_abund = arg
        elif opt in ("-n", "--max_abund"):
            max_abund = arg
        elif opt in ("-b", "--base_range"):
            base_range = arg
            print "Custom base range:", base_range
        elif opt in ("-c", "--chro_list"):
            chro_list = arg
            chro_list = chro_list.split(',')
    
    print '\n\n\n'
    
    logfile(ifile)
            
def hash_smallRNA(chro):
    hash_coords = {}
        
    for line in open(ifile, 'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            token = line.split('\t')
            tokens = token[3].split(';')[1]
            size = int(token[1])
            abund = float(token[2])
            if abund > min_abund and abund < max_abund :
                if base_range == -1:
                    start = 0
                    end = size
                else:
                    start = int(base_range.split('-')[0]) -1
                    end = int(base_range.split('-')[1]) -1
                for i in range(start,end):
                    for coords in tokens.split(','):
                        if len(coords.split('_')) == 2:
                            if len(coords.split('_')[1]) > 0:
                                chromosome = coords.split('_')[0]
                                pos = int(coords.split('_')[1])
                                if chro == chromosome:
                                    if pos+i not in hash_coords:
                                        hash_coords[pos+i] = abund
                                    else:
                                        hash_coords[pos+i] += abund
    return hash_coords

def printOut(chro, hash_coords, o):
    for line in open(dfile,'r'):
        line = line.strip()
        if len(line)>0 and not line.startswith('#'):
            token = line.split('\t')
            chromosome = token[0]
            if chro == chromosome:
                pos = int(token[1])
                o.write(line)
                sum_site = 0
                for i in range(-1*flanking, flanking+1):
                    if i + pos in hash_coords:
                        sum_site += hash_coords[i+pos]
                    
                        
                for i in range(-1*flanking, flanking+1):
                    if i + pos in hash_coords:
                        o.write('\t'+ str(round(hash_coords[i+pos]/float(sum_site),3)))
                    else:
                        o.write('\t'+ '0')
                o.write('\n')
        

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    
    print 'Hold On hashing the small RNA mappings'
    PrinteclapsedTime()
    
    o = open(dfile+'.sRNAmap','w')
    
    for chro in chro_list:
        ### hash the small RNA mappings
        hash_coords = hash_smallRNA(chro)
        
        print chro
        print 'Finished hashing the small RNA mappings'
        print 'Going to print the Site wide tables'
    
        PrinteclapsedTime()
        
        printOut(chro, hash_coords, o)    
        ### close the logfile
    o.close()