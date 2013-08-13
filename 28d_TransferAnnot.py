#-----------------------------------------------------------+
#                                                           |
# 28d_TransferAnnot.py - Transferring GO annotations        |
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
# python ~/script/python/28d_TransferAnnot.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables
global blast, ref, go, HEADER

HEADER = 'database	q.id	q.length	q.frame	q.start	q.end	s.length	s.start	s.end	s.coverage	hsp.coverage	identity	evalue	bit.score	s.title	tax_id	GeneID	status	RNA_nucleotide_accession.version	RNA_nucleotide_gi	protein_accession.version	protein_gi	genomic_nucleotide_accession.version	genomic_nucleotide_gi	start_position_on_the_genomic_accession	end_position_on_the_genomic_accession	orientation	assembly	mature_peptide_accession.version	mature_peptide_gi	Symbol	tax_id	GeneID	GO_ID	Evidence	Qualifier	GO_term	PubMed	Category'

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "28d_TransferAnnot.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 28d_TransferAnnot.py -i <ifile>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global blast, ref, go
    
    blast = ''
    ref = ''
    go = ''
    
    try:
        opts, args = getopt.getopt(argv,"hi:r:g:",["ifile=","ref=","go="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            blast = arg
        elif opt in ("-r", "--ref"):
            ref = arg
        elif opt in ("-g", "--go"):
            go = arg
            
    logfile(blast)
            

def HashFile(file, col_key):            
    hash_file = {}
    for line in open(file, 'r'):
        line = line.strip()
        token = line.split('\t')
        if len(line) > 1 and not line.startswith('#'):
            key = token[col_key -1]
            if key in hash_file:
                hash_file[key] += ';'+ line 
            else:
                hash_file[key] = line
                
    return hash_file


def process_blast(blast, hash_ref, hash_gid):
    print HEADER
    for line in open(blast, 'r'):
        line = line.strip()
        token = line.split('\t')
        
        if len(token[14].split('|')) > 1:
            ref_id = token[14].split('|')[3]
            if ref_id in hash_ref:
                ref_token = hash_ref[ref_id].split(';')
                for i in range(len(ref_token)):
                    if len(ref_token[i].split('\t')) > 1:
                        gene_ID = ref_token[i].split('\t')[1]
                        if gene_ID in hash_gid:
                            tid_token = hash_gid[gene_ID].split(';')
                            for j in range(len(tid_token)):
                                print line + '\t'+ ref_token[i] + '\t' + tid_token[j] 
            
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### hash the gene2go file
    hash_gid = HashFile(go, 2)
   
   
    ### hash the gene2refseq
    hash_ref = HashFile(ref, 6)
    
    ### process the blast results
    process_blast(blast, hash_ref, hash_gid)
    
    ### close the logfile
    o.close()