#-----------------------------------------------------------+
#                                                           |
# 16c_make_hit_table.py                                     |
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
import os,sys,getopt, re


### global variables
global ifile, bblast, pblast

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
    global ifile, bblast, pblast
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:b:p:",["ifile=","blast="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-p", "--pblast"):
            pblast = arg
        elif opt in ("-b", "--bblast"):
            bblast = arg
            
    logfile(ifile)
            
def hash_blast_plant():
    blasthash = {}
    for line in open(pblast,'r'):
        line = line.strip()
        tokens = line.split('\t')
        query = tokens[0]
        evalue = float(tokens[11])
        
        ### header manipulations
        match = re.search(r'_gi\|.+?\|',query)
        if match:
            query = match.group()[1:-1]
        if query.startswith('0'):
            query=query[3:]
        
        if len(query) > 25:
            query = query[0:25]
            
        if query in blasthash:
            if evalue > blasthash[query]:
                blasthash[query] = evalue
        else:
            blasthash[query] = evalue
    
    return blasthash

def hash_blast_bacteria():
    blasthash = {}
    for line in open(bblast,'r'):
        line = line.strip()
        tokens = line.split('\t')
        query = tokens[0]
        evalue = float(tokens[11])
        
        ### header manipulations
        match = re.search(r'_gi\|.+?\|',query)
        if match:
            query = match.group()[1:-1]
        if query.startswith('0'):
            query=query[3:]
        
        if len(query) > 25:
            query = query[0:25]
            
        if query in blasthash:
            if evalue > blasthash[query]:
                blasthash[query] = evalue
        else:
            blasthash[query] = evalue
    
    return blasthash
                            
def parse_proteins(pblasthash, bblasthash):
    HEADER='Accessions\tSymbiosome\tSoybeanNodule\tMedicagoNodule\tLotusRoots\tOrigin\tSpecie'
    print HEADER
    first_line = True
    for line in open(ifile, 'r'):
        if first_line == False:
            line = line.strip()
            tokens = line.split('\t')
            key = tokens[0]
            p_evalue = 0
            b_evalue = -1
            
            if key in pblasthash:
                p_evalue = pblasthash[key]
            if key in bblasthash:
                b_evalue = bblasthash[key]
            
            if p_evalue >= b_evalue:
                print line+'\t'+'P'
            if p_evalue < b_evalue:
                print line+'\t'+'B'
            
                       
        first_line = False

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### hash the blast results
    pblasthash = hash_blast_plant()
    bblasthash = hash_blast_bacteria()
    
    ### parse all the proteins
    parse_proteins(pblasthash, bblasthash)
    
    ### close the logfile
    o.close()