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
    #o = open(ifile+'.burttii_filtered','w')
    first_line = True
    global HEADER
    count = 0
    #print 'Marksers used: ',count
    HEADER='CHROM\tPOS\t'
    string='Burtii_20130605	Gifu_20130609	MG20_genomic_20130609	mg004	mg010	mg012	mg019	\
    mg023	mg036	mg049	mg051	mg062	mg072	mg073	mg077	mg080	mg082	mg083	mg086	mg089\
	    mg093	mg095	mg097	mg101	mg107	mg109	mg112	mg113	mg118	mg123	mg128'
    
    
    for i in string.split():
	HEADER += i+'_GC\t' + i+'_GQ\t' + i+'_DP\t'

    print HEADER

    for line in open(ifile,'r'):
	if len(line) > 1 and not line.startswith('#'):
	    obj = classVCF.VCF(line)
	    genotypes = obj.genotypes()
	    ### check if the MG20 is 0/0 reference Homozygous
	    if obj.genotype(2) != '0/1' and len(obj.alts())==1 and len(obj.refs())==1:
		#o.write(line+'\n')
		count += 1
		genotypes = obj.genotypes()
		
		string = obj.chroms()+'\t'+str(obj.poss()) + '\t'
		for i in range(len(genotypes)):
		    if obj.genotype(i) != 'NONE':
			string += obj.genotype(i) + '\t' + str(obj.genotypeQual(i)) + '\t' +  str(obj.genotypeDepth(i)) + '\t'
		    else:
			string += obj.genotype(i) + '\t' + '0' + '\t' +  '0' + '\t'
		print string
			
 
if __name__ == "__main__":
    
    options(sys.argv[1:])
    empty_file(ifile)
    
    parseFile(ifile)
    
    
    ### close the logfile
    #o.close()