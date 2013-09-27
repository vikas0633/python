#-----------------------------------------------------------+
#                                                           |
# 30a_count_5prime_stacks.py - script for counting 5'       |
# degradome mappings from BAM file                          |
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
# python ~/script/python/30a_count_5prime_stacks.py -i sample.mg010.bam_sorted.bam -r ~/LotusGenome/ljr30/lj_r30.fa --stack_frac_min 0 -c


### import modules
import os,sys,getopt, re
import classCigar


### global variables
global BAMs, ref, MinTotatCov, MinStackCov, StackFraction, HEADER, Genone_coverage

HEADER = 'CHRO\tPOS\tCoverage\tStackSize\tStackFraction'

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "30a_count_5prime_stacks.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 30a_count_5prime_stacks.py
                                            -i <BAMs>
                                            -r <ref>
                                            -t <total_cov_min>
                                            -s <stack_cov_min>
                                            -f <stack_frac_min> 
                                            -c <genone_coverage> [Make mendatory ref.cov file, default: False ]
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global BAMs, ref, MinTotatCov, MinStackCov, StackFraction, Genone_coverage
    BAMs = ''
    ref = ''
    MinTotatCov = 5
    MinStackCov = 5
    StackFraction = 0.5
    Genone_coverage = False
    
    try:
        opts, args = getopt.getopt(argv,"hi:r:t:s:f:c",["BAMs=","ref=","total_cov_min=","stack_cov_min=","stack_frac_min=","genone_coverage="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            BAMs = arg
        elif opt in ("-r", "--ref"):
            ref = arg
        elif opt in ("-t", "--total_cov_min"):
            MinTotatCov = int(arg)
        elif opt in ("-s", "--stack_cov_min"):
            MinStackCov = int(arg)
        elif opt in ("-f", "--stack_frac_min"):
            StackFraction = float(arg)
        elif opt in ("-c", "--genone_coverage"):
            Genone_coverage = True
        
            
class SAM():
    def __init__(self, line):
        tokens = line.split('\t')
        self.id = tokens[0]
        self.flag = tokens[1]
        self.st_chr = tokens[2]
        self.st_pos = tokens[3]
        self.cigar = tokens[5]
        self.read = tokens[9]
            
    
    def __str__(self):
        return self.id
    
    def ids(self):
        return self.id
    
    def flags(self):
        return self.flag
    
    def st_chrs(self):
        return self.st_chr
    
    def st_poss(self):
        return self.st_pos
    
    def reads(self):
        return self.read
    
    def cigars(self):
        return self.cigar
            
def files():
    files = []
    for file in BAMs.split(','):
        print '###', file
        files.append(file.strip())
    return files                        
        
def MakeStrandSpecificSams(file):
    os.system('samtools view -F 0x10 '+ file + ' >forword.temp.sam')
    os.system('samtools view -f 0x10 '+ file + ' >reverse.temp.sam')
    
def GenomeCov(file):
    if not os.path.exists(file+'.cov'):
        os.system('bedtools genomecov -d -ibam '+ file + ' -g '+ ref + ' > '+file+'.cov')
    
def CountPos():
    hash_pos_fw = {}
    hash_pos_rv = {}
    #http://samtools.sourceforge.net/SAMv1.pdf
    # samtools gives the position of leftmost base on genome 
    for line in open('forword.temp.sam','r'):
        line = line.strip()
        obj = SAM(line)
        key = obj.st_chrs() , obj.st_poss()
        if key not in hash_pos_fw:
            hash_pos_fw[ key ] = 1
        else:
            hash_pos_fw[ key ] += 1
        
    ### for minus strand 5' end of the read is the pos + len(read) -1
    for line in open('reverse.temp.sam','r'):
        line = line.strip()
        obj = SAM(line)
        cigar = classCigar.Cigar(obj.cigars())
        ins = int(cigar.insertions())
        dele = int(cigar.deletions())
        
        key = obj.st_chrs() , str(int(obj.st_poss()) + len(obj.reads()) - 1 + dele - ins)
        if key not in hash_pos_rv:
            hash_pos_rv[ key ] = 1
        else:
            hash_pos_rv[ key ] += 1
    
    
    return hash_pos_fw, hash_pos_rv

def print_out(file, var, out):
    ### loop though the genome refenece count
    for line in open(file, 'r'):
        line = line.strip()
        tokens = line.split('\t')
        cov = int(tokens[2])
        key = tokens[0],tokens[1]
        if key in var:
            if int(cov) > MinTotatCov:# check total coverage
                if int(var[key]) > MinStackCov: # check stack coverage
                    stack_fraction = var[key]/float(cov)
                    if stack_fraction > StackFraction: # check stack fraction
                        
                        out.write(str(tokens[0])+'\t'+str(tokens[1])+'\t'+str(cov) +'\t'+ str(var[key])+'\t'+ str(round(stack_fraction,2))+'\n')
    
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### get the file list
    file_list = files()
    
    ### make forward and reverse read containing SAM files
    for file in file_list:
        
        ### make a genome coverage file
        if Genone_coverage == True:
            GenomeCov(file)
        
        MakeStrandSpecificSams(file)    
        out = open(file+'.sites','w')
        out.write(HEADER+'\n')
        ### count the positions on the genome
        hash_pos_fw, hash_pos_rv = CountPos()
        
        ### print the position passing criteria
        print_out(file+".cov", hash_pos_fw, out)
        print_out(file+".cov", hash_pos_rv, out)
        out.close()
    ### close the logfile
    o.close()