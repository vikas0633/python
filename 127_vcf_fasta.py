#-----------------------------------------------------------+
#                                                           |
# 127_vcf_fasta.py - Script to convert VCF to fasta file    |
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
global ifile

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

def printFasta(seqs, genotypes_names):
    o = open(ifile+'.fasta','w')
    for i in range(len(seqs)):
        o.write('>'+genotypes_names[i]+'\n')
        o.write(seqs[i]+'\n')
    
    o.close()


def makeFasta():
    seqs = {}
    count = 0
    for line in open(ifile,'r'):
        if line.startswith('#CHROM'):
            line = line.strip('\n')
            obj = classVCF.VCF(line)
            genotypes_names = obj.genotypes()
            
        
        
        if len(line) > 1 and not line.startswith('#'):
            line = line.strip('\n')
            obj = classVCF.VCF(line)
            genotypes = obj.genotypes()
            
            count += 1
            if count%100000 == 0:
                print 'Number of markers processed: ', '{:9,.0f}'.format(count)

            ### check if the burttii is not 0/1 heterogyzous
            for i in range(len(genotypes)):
                if i in seqs:
                    if obj.genotype(i) == '0/0':
                        seqs[i] += obj.refs()
                    elif obj.genotype(i) == '1/1':
                        seqs[i] += obj.alts()
                    else:
                        seqs[i] += 'N'
                else:
                    if obj.genotype(i) == '0/0':
                        seqs[i] = obj.refs()
                    elif obj.genotype(i) == '1/1':
                        seqs[i] = obj.alts()
                    else:
                        seqs[i] = 'N'
                            
    printFasta(seqs, genotypes_names)

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    makeFasta()
    
    ### close the logfile
    o.close()