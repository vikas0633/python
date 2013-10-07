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
import os,sys,getopt, re, classVCF


### global variables
global ifile, window_size, step_size, window_count, step_count, callable_file, count_based

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
                                        -j <callable_file>
                                        -w <window_size_pos> ### Position based
                                        -c <window_size_count> ### count based
                                        
                                        
            File must follow the format:
            <chromosome> <Position> <DataPoints>
            
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile, window_size, step_size, window_count, step_count, callable_file, count_based
    ifile = ''
    callable_file = ''
    window_size = 1000
    step_size = 100
    window_count = 1000
    step_count = 100
    count_based = False
    try:
        opts, args = getopt.getopt(argv,"hi:j:w:c:s:d:a",["ifile=","callable_file=","window_size_pos=","step_size_pos=","window_size_count=","step_size_count=","alt" ])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-j", "--callable_file"):
            callable_file = arg
        elif opt in ("-w", "--window_size_pos"):
            window_size = int(arg)
        elif opt in ("-c", "--window_size_count"):
            window_count = int(arg)
        elif opt in ("-s", "--step_size_pos"):
            step_size = int(arg)
        elif opt in ("-d", "--step_size_count"):
            step_count = float(arg)
        elif opt in ("-a", "--alt"):
            count_based = True
        
    logfile(ifile)
            
def get_size(file):
    count = 0
    chroHash = {}
    for line in open(file,'r'):
        count += 1
        ### print the lines processed
        if count%100000 == 0:
            print 'Lines processed: ', '{:9,.0f}'.format(count)
        line = line.strip()
        if len(line)>1 and not line.startswith('#'):
            token = line.split('\t')
            if token[0] not in chroHash:
                chroHash[token[0]] = '' 				
    return chroHash

def parse_call(chromosome,callable_file):
    hash = {}
    count = 0
    for line in open(callable_file,'r'):
        if len(line)>1 and not line.startswith('#'):
            #obj = classVCF.VCF(line)
            count += 1
            ### print the lines processed
            if count%100000 == 0:
                print 'Lines processed in callable fraction file: ', chromosome, '{:9,.0f}'.format(count)
            line = line.strip()
            tokens = line.split('\t')
            if tokens[0]==chromosome:
                for pos in range(int(tokens[1]), int(tokens[2])+1):
                    hash[pos] = ''
    return hash
 
def parse(chromosome, o):
    count = 0
    data_hash = {}
    for line in open(ifile, 'r'):
        if len(line)>1 and not line.startswith('#'):
            count += 1
            ### print the lines processed
            if count%100000 == 0:
                print 'Lines processed in the variant file: ', chromosome, '{:9,.0f}'.format(count)
            line = line.strip()
            tokens = line.split('\t')
            if tokens[0]==chromosome:
                pos = int(tokens[1])
                data_hash[pos] = float(tokens[2]), float(tokens[3]), float(tokens[4]), float(tokens[5]), float(tokens[6]), float(tokens[7]) ,\
                float(tokens[8]), float(tokens[9])
                
                
    return data_hash
def window(chromosome, o, data_hash, call_hash):
    ### run the window position based analysis    

    if len(data_hash)>0:
        start = min(data_hash)
        end = max(data_hash)

    for i in range(start, end+1, step_size):
        print 'Lines processed printing out: ', chromosome, '{:9,.0f}'.format(i)
        pos_sum = 0
        snp_sum = 0
        depth_sum = 0
        genotypeCalls_sum = 0
        genotypeCallsHete_sum = 0
        genotypeCallsHomo_sum = 0
        InbreedingCoeffs_sum = 0
        HaplotypeScores_sum = 0
        temp_count = 0
        for j in range(i, i+window_size):
            pos_sum += 1
            if j in call_hash:
                temp_count += 1
                depth_sum += 0
            if j in data_hash:
                call, snp, depth, genotypeCalls, genotypeCallsHete, genotypeCallsHomo, InbreedingCoeffs, HaplotypeScores = data_hash[j]
                snp_sum += snp
                genotypeCalls_sum += genotypeCalls
                genotypeCallsHete_sum += genotypeCallsHete
                genotypeCallsHomo_sum += genotypeCallsHomo
                InbreedingCoeffs_sum += InbreedingCoeffs
                HaplotypeScores_sum += HaplotypeScores
        
        
            
        if temp_count > 0:
        
            if snp_sum > 0:
                o.write(chromosome+'\t'+ \
                        str(i+window_size/2) + '\t' + \
                str(temp_count/float(pos_sum))+ '\t' + \
                str(snp_sum/float(temp_count))+ '\t' + \
                str(depth_sum/float(temp_count))+ '\t'+ \
                str(genotypeCalls_sum/float(snp_sum))+ '\t' + \
                str(genotypeCallsHete_sum/float(snp_sum))+ '\t'+ \
                str(genotypeCallsHomo_sum/float(snp_sum))+ '\t'+ \
                str(InbreedingCoeffs_sum/float(snp_sum))+ '\t'+ \
                str(HaplotypeScores_sum/float(snp_sum))+ '\n')
            else:
                o.write(chromosome+'\t'+ \
                        str(i+window_size/2) + '\t' + \
                str(temp_count/float(pos_sum))+ '\t' + \
                str(snp_sum/float(temp_count))+ '\t' + \
                str(depth_sum/float(temp_count))+ '\t'+ \
                str(0)+ '\t' + \
                str(0)+ '\t'+ \
                str(0)+ '\t'+ \
                str(0)+ '\t'+ \
                str(0)+ '\n')

def window_call(chromosome, o, data_hash, call_hash):
    if len(data_hash)>0:
        start = min(data_hash)
        end = max(data_hash)
    new_block = True
    pos_sum = 0
    snp_sum = 0
    depth_sum = 0
    genotypeCalls_sum = 0
    genotypeCallsHete_sum = 0
    genotypeCallsHomo_sum = 0
    InbreedingCoeffs_sum = 0
    HaplotypeScores_sum = 0
    temp_count = 0
    w_count = 0
    for i in range(start, end+1):
        pos_sum += 1
        if i%1000000 == 0:
            print 'Lines processed printing out: ', chromosome, '{:9,.0f}'.format(i)
            
        if i in call_hash:
            w_count += 1
        
        if i in data_hash:
            call, snp, depth, genotypeCalls, genotypeCallsHete, genotypeCallsHomo, InbreedingCoeffs, HaplotypeScores = data_hash[i]
            snp_sum += snp
            genotypeCalls_sum += genotypeCalls
            genotypeCallsHete_sum += genotypeCallsHete
            genotypeCallsHomo_sum += genotypeCallsHomo
            InbreedingCoeffs_sum += InbreedingCoeffs
            HaplotypeScores_sum += HaplotypeScores
        
        if w_count == window_count:
            new_block = True
            if snp_sum > 0:
                o.write(chromosome+'\t'+ \
                        str((start+i)/2) + '\t' + \
                str(w_count/float(pos_sum))+ '\t' + \
                str(snp_sum/float(w_count))+ '\t' + \
                str(depth_sum/float(w_count))+ '\t'+ \
                str(genotypeCalls_sum/float(snp_sum))+ '\t' + \
                str(genotypeCallsHete_sum/float(snp_sum))+ '\t'+ \
                str(genotypeCallsHomo_sum/float(snp_sum))+ '\t'+ \
                str(InbreedingCoeffs_sum/float(snp_sum))+ '\t'+ \
                str(HaplotypeScores_sum/float(snp_sum))+ '\n')
            
        else:
            new_block = False
        
        if new_block == True:
            start = i
            w_count = 0
            pos_sum = 0
            snp_sum = 0
            depth_sum = 0
            genotypeCalls_sum = 0
            genotypeCallsHete_sum = 0
            genotypeCallsHomo_sum = 0
            InbreedingCoeffs_sum = 0
            HaplotypeScores_sum = 0
        
        

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    HEADER = '#Chromosome\tPosition\tCallable\tSNP\tDepth\tgenotypeCalls\tgenotypeCallsHete\tgenotypeCallsHomo\tInbreedingCoeffs\tHaplotypeScores'
    
    print 'Hashing the chromosomes name'
    chroHash = get_size(ifile)
    
    o = open(ifile+'.w'+str(window_size)+'.s'+str(step_size), 'w')
    o.write(HEADER+'\n')
    
    for chromosome in sorted(chroHash):
        if chromosome != 'chr0':
            ### hash the callable fraction
            call_hash = parse_call(chromosome,callable_file)
            
            ### parse file
            data_hash = parse(chromosome, o)
            if count_based == False:
                window(chromosome, o, data_hash, call_hash)
            else:
                window_call(chromosome, o, data_hash, call_hash)
    
    ### close the logfile
    o.close()