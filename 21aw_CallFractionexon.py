#-----------------------------------------------------------+
#                                                           |
# 21aw_CallFractionexon.py -                                |
# script to calculate the callable fraction by transcripts  |
# in exonic regions                                         |
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
import os,sys,getopt, re, E_get_chr_size_gff3, classGene, classmRNA


### global variables
global ifile, call

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
            python 100b_fasta2flat.py -i <ifile> ### GFF3 file
                                        -c <callable> ### callable fraction coords
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile, call
    ifile = ''
    call = ''
    try:
        opts, args = getopt.getopt(argv,"hi:c:",["ifile=","call="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-c", "--call"):
            call = arg
            
    logfile(ifile)
            
def parseCall(chrom):
    hash_call = {}
    for line in open(call, 'r'):
        if len(line) > 1 and not line.startswith('#'):
            line = line.strip()
            token = line.split()
            hash_call[int(token[1])] = ''
        
    return hash_call
                    
                            
def get_exon_fraction(chrom, hash_call):
    HEADER = 'Lj30_ID\tExonLength\tCDSLength\tCallableExon\tCallableCDS'
    obj_list=[]
    for line in open(ifile, 'r'):
        if len(line) > 1 and not line.startswith('#'):
            line = line.strip()
            obj = classGene.GFF3(line)
            if obj.types() == "mRNA":
                obj_mRNA = classmRNA.mRNA(line, obj)
                obj_list.append(obj_mRNA)
            if obj.types() == "mRNA" or obj.types() == "exon" or obj.types() == "CDS":
                obj_mRNA.AddData(line, obj)
    print HEADER
    for obj_mRNA in obj_list:
        print str(obj_mRNA), obj_mRNA.GetExonLength(), obj_mRNA.GetCDSLength(), obj_mRNA.GetExonicOverlap(hash_call), obj_mRNA.GetCDSOverlap(hash_call)

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    size = E_get_chr_size_gff3.get_size(ifile)
    
    for chrom in size:
        hash_call = parseCall(chrom)
        get_exon_fraction(chrom, hash_call)
    
    ### close the logfile
    o.close()