#-----------------------------------------------------------+
#                                                           |
# 112_iprscanout_parser.py - IPRScan XML parser             |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# IPRScan XML parser                                        |
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


'''
def IPR:
    
    def __ini__(self, block):
'''        

def header():
    print 'Gene_ID' +'\t'+ \
    'Transcript_ID' + '\t' + \
    'protein_ID' + '\t' + \
    'GO_ID' + '\t' + \
    'Category' + '\t' + \
    'Description'
     


def get_IDs(line):
    match = re.search(r'\"(.+?)\"',line)
    if match:
        protein_ID = match.group().replace('"','')
        transcript_id = '_'.join(protein_ID.split('_')[:-2]).replace('"','')
        gene_id = '.'.join(transcript_id.split('.')[:-1]).replace('"','')
    
    return gene_id, transcript_id, protein_ID

def get_GOid(line):
    match = re.search(r'\"(.+?)\"',line)
    if match:
        GO_id = match.group().replace('"','')
    return GO_id

def parse(infile):
    
    for line in open(infile,'r'):
        line = line.strip()
    
        if line.startswith("<protein id"):
            g_id, t_id, p_id = get_IDs(line)
                
        if line.startswith("<classification"):
            GO_id = get_GOid(line)
        if line.startswith("<category>"):
            cat = line.replace('<category>', '').replace('</category>','')
        if line.startswith('<description>'):
            desc = line.replace('<description>','').replace('</description>','')
        if line.startswith('</classification>'):
            if GO_id !='':
                print g_id +'\t'+ t_id +'\t'+ p_id +'\t'+ GO_id +'\t'+ cat +'\t'+ desc 
    
### main argument to 
def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
            
    return infile
    
                            
        

if __name__ == "__main__":
    
    infile = options(sys.argv[1:])
    
    ### print header
    header()
    
    ###parse the IPRscan output
    parse(infile)
    
    ### close the logfile
    o.close()