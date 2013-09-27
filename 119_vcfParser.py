#-----------------------------------------------------------+
#                                                           |
# 119_vcfParser.py - script to parse vcf format file       |
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
# python ~/Desktop/script/python/119_vcfParser.py -i snp.90.PhredQual_5000.vcf


### import modules
import os,sys,getopt, re, classVCF


### global variables
global ifile, HEADER

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
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
            
    logfile(ifile)
            
### check if file empty
def empty_file(infile):
    if os.stat(infile).st_size==0:
        sys.exit('File is empty')
                                
        
def parseFile(ifile):
    o = open(ifile+'.MG20filtered','w')
    global HEADER
    count = 0
    for line in open(ifile,'r'):
        if len(line) > 1 and not line.startswith('##'):
            line = line.strip('\n')
            if line.startswith('#CHROM'):
                o.write(line+'\n')
                HEADER = line
                samples_het = []
                samples_homo = []
                sample_names = line.split('\t')[9:]
                samples_len = len(line.split('\t')) -9
                for i in range(samples_len):
                    samples_het.append(0)
                    samples_homo.append(0)
            else:
                obj = classVCF.VCF(line)
                genotypes = obj.genotypes()
                ### check if the MG20 is 0/0 reference Homozygous
                if obj.genotype(2) == '0/0':
                    o.write(line+'\n')
                    count += 1
                    genotypes = obj.genotypes()
                    for i in range(len(genotypes)):
                        if obj.genotype(i) =='0/1' or obj.genotype(i) =='1/0':
                            samples_het[i] += 1
                        elif obj.genotype(i) =='0/0' or obj.genotype(i) =='1/1':
                            samples_homo[i] += 1
    print 'Marksers used: ',count
    print 'Sample\tHetCount\tHomoCount\tHetPer\tHomoPer'
    
    for i in range(len(sample_names)):
        total = int(samples_het[i]) + int(samples_homo[i])
        Het_per = float(samples_het[i])/total
        Homo_per = float(samples_homo[i])/total
        print sample_names[i] + '\t' + str(samples_het[i]) + '\t' + str(samples_homo[i]) + '\t' + str(Het_per) + '\t' + str(Homo_per)
    o.close()
 
if __name__ == "__main__":
    
    options(sys.argv[1:])
    empty_file(ifile)
    
    parseFile(ifile)
    
    
    ### close the logfile
    o.close()