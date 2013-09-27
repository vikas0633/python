#-----------------------------------------------------------+
#                                                           |
# 115_MapFastq.py - Script to Map Fastq files               |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/115_MapFastq.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os, sys, getopt, re, glob


### global variables
global ref, folder, extension, fastqc, mapper, threads, seed_length, mismatches, single, index, compress_extension, uncompress

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "115_MapFastq.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 115_MapFastq.py
                                    -r <ref> [reference sequence]
                                    -f <folder> [folder1, folder2, .., folderN]
                                    -x <extension> [default: fastq]
                                    -q <fastqc> [#runs fastq rather than mapping]
                                    -p <mapper> [default: bwa]
                                    -t <threads> [default: 6, numbers of core to be used]
                                    -l <seed_length> [default: 28, seed length to be used in mapping]
                                    -m <mismatches> [default: 2, mismatches allowed in the seed]
                                    -s <single> [default: Pair End alignments]
                                    -i <index> [Option to create index for Reference Sequence]
                                    -u <uncompress> [default extention: bz2]
                                    
            Fastq pair must be specified with "*_R1_*" and "*_R2_*"  
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ref, folder, extension, fastqc, mapper, threads, seed_length, mismatches, single, index, compress_extension, uncompress
    ref = ''
    folder = ''
    extension = 'fastq'
    fastqc = False
    mapper='bwa'
    threads = 6
    seed_length = 28
    mismatches = 2
    single = False
    index = False
    compress_extension = 'bz2'
    uncompress = False
    
    try:
        opts, args = getopt.getopt(argv,"hr:f:x:qp:t:l:m:s:iu:",["ref=","folder=","extension=","fastqc=","mapper=","threads=","seed_length=","mismatches=","single=","index=","uncompress="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-r", "--ref"):
            ref = arg
        elif opt in ("-f", "--folder"):
            folder = arg
        elif opt in ("-x", "--extension"):
            extension = arg
        elif opt in ("-q", "--fastqc"):
            fastqc = True
        elif opt in ("-p", "--mapper"):
            mapper = arg
        elif opt in ("-t", "--threads"):
            threads = arg
        elif opt in ("-l", "--seed_length"):
            seed_length = arg
        elif opt in ("-m", "--mismatches"):
            mismatches = arg
        elif opt in ("-s", "--single"):
            single = True
        elif opt in ("-i", "--index"):
            index = True
        elif opt in ("-u", "--uncompress"):
            uncompress = True
            compress_extension = arg
    
    logfile(ref)
            
    return  
    
def Uncompress(file):
    
    if compress_extension == 'bz2':
        os.system('bzip2 -d --keep --verbose ' + file)
    
def files():
    print 'Files to processed'
    file_list = []
    for f in folder.split(','):
        f = f.strip()
        if uncompress == True:
            for file in glob.glob(os.path.join(f, '*'+compress_extension)):
                Uncompress(file) 
                file_list.append('.'.join(file.split('.')[:-1]))
        else:
            for file in glob.glob(os.path.join(f, '*'+extension)):
                file_list.append(file)
                
        print file
        
    return file_list

def FastQC(file_list):
    if fastqc == True:
        for file in file_list:
            os.system('Running FastQC for '+file)
            os.system('fastqc '+file)

def Index():
    if index == True:
        os.system('nice -n 19 samtools faidx '+ ref)
        os.system('nice -n 19 bwa index -a bwtsw '+ ref)
    if not os.path.isfile(ref + '.fai'):
        os.system('nice -n 19 samtools faidx '+ ref)
        
        
def AlignReads(file_list):
    for file in file_list:
        if mapper == 'bwa':
            if not os.path.isfile(file+'.sai'):
                os.system('nice -n 19 bwa aln -t '+ str(threads) +' -l ' +str(seed_length)+ ' ' + ref + ' ' + file + ' > ' + file+'.sai')
    
def MapReads(file_list):
    for file in file_list:
        if mapper == 'bwa':
            if single == False:
                if re.search('_R1_', file):
                    read1 = file
                    read2 = file.replace('_R1_','_R2_')
                    rg = file.strip().split('/')[-1].strip()[:6].strip()
                    rg = '"@RG\tID:'+rg+'\tSM:'+rg+'\tPL:illumina\tLB:lib1\tPU:unit"'
                    
                    os.system('nice -n 19 bwa sampe -P '+ ' -r ' + rg +' '+ ref +' '+ read1+".sai " + read2+".sai " +\
                              read1 +' '+read2 +' | nice -n 19 samtools view -bt '+ ref+'.fai -| nice -n 19 samtools sort - '+ \
                              read1+'_sorted')
                    
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    
    ### print the files to be process and return the file path as list
    file_list = files()
    
    ### fastqc
    FastQC(file_list)
    
    ### index reference
    Index()
    
    ### align the reads
    AlignReads(file_list)
    
    ### map the files
    MapReads(file_list)
    
    ### close the logfile
    o.close()