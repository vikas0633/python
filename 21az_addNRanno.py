#-----------------------------------------------------------+
#                                                           |
# 21az_addNRanno.py - script to add nr annotation to GFF3 file  |
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
import os,sys,getopt, re, classGene


### global variables
global ifile, anno

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
            python 100b_fasta2flat.py -i <ifile> -a <anno>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile, anno
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:a:",["ifile=","anno="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-a", "--anno"):
            anno = arg
            
    logfile(ifile)
            
def HashAnno():
    hash_anno = {}
    
    for line in open(anno,'r'):
        line = line.strip()
        if len(line)>0 and not line.startswith('#'):
            tokens = line.split('\t')
            hash_anno[tokens[0]] = tokens[1]
        
    return hash_anno
                            
        
def parse(hash_anno):
    
    for line in open(ifile,'r'):
        line = line.strip()
        if len(line)>0 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            
            if obj.types() == 'mRNA':
                print line + ';'+'Annotation="'+hash_anno[str(obj)].replace(',','')+'"'
            else:
                print line




if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### hash anno
    hash_anno = HashAnno()
    
    ### GFF3 file
    parse(hash_anno)
    
    ### close the logfile
    o.close()