#-----------------------------------------------------------+
#                                                           |
# 30c_degradome_site_miRNAs.py - script for counting miRNA spicies on cleavage site      |
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
global ifile, time_start, dfile, chro_list
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
            python 100b_fasta2flat.py -i <ifile>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile, dfile, chro_list
    ifile = ''
    dfile = ''
    chro_list = ['chr0','chr1','chr2','chr3','chr4','chr5','chr6','chloro','mito']
    try:
        opts, args = getopt.getopt(argv,"hi:a:",["ifile=","dfile="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-a", "--dfile"):
            dfile = arg
        elif opt in ("-c", "--chro_list"):
            chro_list = arg
            chro_list = chro_list.split(',')
    logfile(ifile)
            
    
def hash_miRNA(chro):
    hash_coords = {}
        
    for line in open(ifile, 'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            token = line.split('\t')
            tokens = token[3].split(';')[1]
            size = int(token[1])
            abund = float(token[2])
            seq = token[0]
            for i in range(0,size):
                for coords in tokens.split(','):
                    if len(coords.split('_')) == 2:
                        if len(coords.split('_')[1]) > 0:
                            chromosome = coords.split('_')[0]
                            pos = int(coords.split('_')[1])
                            if chro == chromosome:
                                if pos+i not in hash_coords:
                                    hash_coords[pos+i] = seq
                                else:
                                    hash_coords[pos+i] += ',' + seq

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
                if pos in hash_coords:
                    o.write('\t'+ str(len(set(hash_coords[pos].split(','))))+'\n')
                else:
                    o.write('\t'+ '0'+'\n')
            


if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    print 'Hold On hashing the miRNA  mappings'
    PrinteclapsedTime()
    
    o = open(dfile+'.miRNAmap','w')
    
    for chro in chro_list:
        ### hash the small RNA mappings
        hash_coords = hash_miRNA(chro)
        
        print chro
        print 'Finished hashing the micro RNA mappings'
        print 'Going to print the Site wide tables'
    
        PrinteclapsedTime()
        
        printOut(chro, hash_coords, o)    
        ### close the logfile
    o.close()