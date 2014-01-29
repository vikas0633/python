#-----------------------------------------------------------+
#                                                           |
# 121_SlidingWindow.py - plot data with sliding window on multiple chromosome |
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
import os,sys,getopt, re, threading
from multiprocessing import Process, Queue

### global variables
global ifile, window_size, step_size, window_count, step_count

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 100b_fasta2flat.py -i <ifile>
                                        -w <window_size_pos> ### Position based
                                        -c <window_size_count> ### count based
                                        
                                        
            File must follow the format:
            <chromosome> <Position> <Group>
            
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile, window_size, step_size, window_count, step_count
    ifile = ''
    window_size = 1000000
    step_size = 100000
    window_count = 1000
    step_count = 100
    try:
        opts, args = getopt.getopt(argv,"hi:w:c:s:d:",["ifile=","window_size_pos=","step_size_pos=","window_size_count=","step_size_count=" ])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-w", "--window_size_pos"):
            window_size = int(arg)
        elif opt in ("-c", "--window_size_count"):
            window_count = float(arg)
        elif opt in ("-s", "--step_size_pos"):
            step_size = int(arg)
        elif opt in ("-d", "--window_size_count"):
            step_count = float(arg)
            
    logfile(ifile)
            
def get_size(file):
    first_line = True
    count = 0
    chroHash = {}
    for line in open(file,'r'):
        count += 1
        ### print the lines processed
        if count%100000 == 0:
            print 'Lines processed: ', '{:9,.0f}'.format(count)
        line = line.strip()
        if len(line)>1 and not line.startswith('#') and first_line == False:
            token = line.split('\t')
            if token[0] not in chroHash:
                chroHash[token[0]] = ''
        first_line = False
    return chroHash

def window(chromosome, o, data_hash):
    ### run the window position based analysis    
    
    
    start = 0
    end = 200000000
    
    '''
    if len(data_hash)>0:
        start = min(data_hash)
        end = max(data_hash)
   '''
    for i in range(start, end+1, step_size):
        print 'Lines processed printing out: ', chromosome, '{:9,.0f}'.format(i)
        temp_sum = 0
        for j in range(i, i+window_size):
            if j in data_hash:
                temp_sum += data_hash[j]
        o.write(chromosome+'\t'+ str(i+window_size/2) + '\t' + str(temp_sum)+'\n')
            
def parse(chromosome, count):
    #global data_hash, data_type
    o = open('tempThread'+str(count), 'w')
    first_line = False
    count = 0
    data_hash = {}
    for line in open(ifile, 'r'):
        if len(line)>1 and not line.startswith('#') and first_line == False:
            count += 1
            ### print the lines processed
            if count%100000 == 0:
                print 'Lines processed: ', chromosome, '{:9,.0f}'.format(count)
            line = line.strip()
            tokens = line.split('\t')
            if tokens[0]==chromosome:
                pos = int(tokens[1])
                if pos in data_hash:
                    data_hash[pos] += 1
                else:
                    data_hash[pos] = 1
        first_line = False
    window(chromosome, o, data_hash)


if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    HEADER = 'Chromosome\tPosition\tData'
    o.write(HEADER+'\n')
    
    print 'Hashing the chromosomes name'
    chroHash = get_size(ifile)
    
  
    ### multithreading
    thread_list = []
    count = 0
    for chromosome in sorted(chroHash):
        count += 1
        t = Process(target=parse, args=(chromosome,count,))
        t.start()
        thread_list.append(t)
    for thread in thread_list:
        thread.join()
        
    '''
    for chromosome in sorted(chroHash):
        ### parse file
        #data_hash,data_type  = parse(chromosome, o)
        window(chromosome, o, data_hash, data_type)
    '''
    ### close the logfile
    os.system('cat tempThread* >' +ifile+'.w'+str(window_size)+'.s'+str(step_size))
    os.system('rm tempThread*')
    