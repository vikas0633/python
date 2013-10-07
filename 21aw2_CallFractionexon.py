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
        if len(line) > 1 and line.startswith(chrom):
            line = line.strip()
            token = line.split()
            hash_call[int(token[1])] = ''
        
    return hash_call
                    
                            
def get_exon_fraction(chrom, hash_call):
    first_transcript = True
    hash_exon = {}
    hash_cds = {}
    HEADER = 'Lj30_ID\tExonLength\tCDSLength\tCallableExon\tCallableCDS'
    print HEADER
    for line in open(ifile, 'r'):
        if len(line) > 1 and not line.startswith('#'):
            line = line.strip()
            obj = classGene.GFF3(line)
            if obj.seqids() == chrom:
                if obj.types() == "mRNA":
                    if first_transcript == False:
                        exon_len = len(hash_exon)
                        cds_len = len(hash_cds)
                        exon_call_len = 0
                        cds_call_len = 0
                        for i in hash_exon:
                            if i in hash_call:
                                exon_call_len += 1
                        for i in hash_cds:
                            if i in hash_call:
                                cds_call_len += 1
                        print id+'\t'+str(exon_len)+'\t'+str(exon_call_len)+'\t'+str(cds_len)+'\t'+str(cds_call_len)
                    first_transcript = False
                    hash_exon = {}
                    hash_cds = {}
                    id = str(obj)
                elif obj.types() == "exon":
                    for i in range(int(obj.starts()), int(obj.ends())):
                        hash_exon[i] = ''
                elif obj.types() == "CDS":
                    for i in range(int(obj.starts()), int(obj.ends())):
                        hash_cds[i] = ''
    print id+'\t'+str(exon_len)+'\t'+str(exon_call_len)+'\t'+str(cds_len)+'\t'+str(cds_call_len)
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    size = E_get_chr_size_gff3.get_size(ifile)
    
    for chrom in size:
        hash_call = parseCall(chrom)
        get_exon_fraction(chrom, hash_call)
    
    ### close the logfile
    o.close()