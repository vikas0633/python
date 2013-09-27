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
global ifile, blast

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
    global ifile, blast
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:b:",["ifile=","blast="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-b", "--blast"):
            blast = arg
            
    logfile(ifile)
            
def hash_blast():
    blasthash = {}
    list_sample=['01','04','05','07']
    for line in open(blast,'r'):
        line = line.strip()
        tokens = line.split('\t')
        query, hit = tokens[0], tokens[1]
        if query[:2] in list_sample and hit[:2] in list_sample:
            if len(query) > 20:
                query = query[0:20]
            
            if query not in blasthash:
                blasthash[query] = [query[:2],hit[:2]]
            else:
                blasthash[query].append(hit[:2])
            
            if len(hit) > 20:
                hit = hit[0:20]
            if hit not in blasthash:
                blasthash[hit] = [query[:2],hit[:2]]
            else:
                blasthash[hit].append(query[:2])
    return blasthash
                            
def parse_proteins(blasthash):
    HEADER='Accessions\tSymbiosome\tSoybeanNodule\tMedicagoNodule\tLotusRoots'
    print HEADER
    hash_data={}
    list_sample=['01','04','05','07']
    for line in open(ifile, 'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            if line.startswith('>'):
                header = line[1:]
                if header[:2] in list_sample:
                    header = header.split(':')[0].strip()
                    if len(header) > 20:
                        header = header[0:20]
                    
                    samples = []
                    
                    if header in blasthash:
                        for sample in list_sample:
                            if sample in blasthash[header]:
                                samples.append('+')
                            else:
                                samples.append('-')
                    else:
                        sam = header[:2]
                        for sample in list_sample:
                            if sample != sam:
                                samples.append('-')
                            else:
                                samples.append('+')
                                
                    match = re.search(r'_gi\|.+?\|',header)
                    if match:
                        header = match.group()[1:-1]
                    
                    string = header.replace('>','')
                    if header.startswith('0'):
                        header=header[3:]
                    
                    if header not in hash_data:
                        hash_data[header] = samples
                    else:
                        for i in range(len(list_sample)):
                            if hash_data[header][i] != samples[i]:
                                hash_data[header][i] = '+'
    for acc in hash_data:
        print acc+'\t'+'\t'.join(hash_data[acc])

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### hash the blast results
    blasthash = hash_blast()
    
    ### parse all the proteins
    parse_proteins(blasthash)
    
    ### close the logfile
    o.close()