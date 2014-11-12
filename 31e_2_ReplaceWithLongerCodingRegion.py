#-----------------------------------------------------------+
#                                                           |
# 31e_2_ReplaceWithLongerCodingRegion.py - Script to find   |
# the longest protein coding evidence with overlapping exons|
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
import classGene

### global variables
global infile

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
            python 31e_2_ReplaceWithLongerCodingRegion.py
                                                        -r <replacement_file>
                                                        -c <combined_gff3>
                                                        -l <last_gff3>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global replacement_file, combined_gff3, last_gff3
    try:
        opts, args = getopt.getopt(argv,"hr:c:l:",["replacement_file=","combined_gff3=","last_gff3="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-r", "--replacement_file"):
            replacement_file = arg
        elif opt in ("-c", "--combined_gff3"):
            combined_gff3 = arg
        elif opt in ("-l", "--last_gff3"):
            last_gff3 = arg

    
def hash_replacementIds():
    '''
    Lj0g0000369.1	model.Ljchr0_pseudomol_20120828.path1.gene55	GlimmerHMM	chr0	147185	156	151
    Lj0g0000509.1	model.Ljchr0_pseudomol_20120828.path1.gene76	GlimmerHMM	chr0	223430	318	0
    '''
    reaplcement_IDs = {}
    hash_combine_ID = {}
    for line in open(replacement_file, 'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            token = line.split('\t')
            reaplcement_IDs[token[0]] = token[1] ### storing corresponding ids
            hash_combine_ID[token[1]] = ''
    
    return reaplcement_IDs, hash_combine_ID

def hash_combined(hash_combine_ID):
    hash_combine = {}
    count = 0
    for line in open(combined_gff3, 'r'):
        line = line.strip()
        count += 1
        if re.search('ID=',line):
            if len(line) > 0 and not line.startswith('#'):
                if count%100000 == 0:
                    print 'Hashing combined: ', '{:9,.0f}'.format(count)
                obj = classGene.GFF3(line)
                if obj.types() == 'mRNA':
                    ID = str(obj)
                    if ID in hash_combine_ID:
                        hash_combine[ID] = {}
                        hash_combine[ID]['mRNA'] = line
                        hash_combine[ID]['exon'] = ''
                        hash_combine[ID]['CDS'] = ''
                        hash_lines = True
                    else:
                        hash_lines = False
                elif obj.types() == 'exon':
                    if hash_lines == True:
                        hash_combine[ID]['exon'] += (','+line)
                elif obj.types() == 'CDS':
                    if hash_lines == True:
                        hash_combine[ID]['CDS'] += (','+line)
    return hash_combine

def make_gff3(reaplcement_IDs, hash_combine):
    out = open(last_gff3+'.replaced','w')
    count = 0
    for line in open(last_gff3, 'r'):
        line = line.strip()
        count += 1
        if len(line) > 0 and not line.startswith('#'):
            if count%100000 == 0:
                print 'Printing final GFF3: ', '{:9,.0f}'.format(count)
            obj = classGene.GFF3(line) 
            if obj.types() == "gene":
                out.write(line+'\n')
                source = obj.sources()
                g_id = str(obj)
            elif  obj.types() == "mRNA":
                ID = str(obj)
                if ID in reaplcement_IDs:
                    token = hash_combine[reaplcement_IDs[ID]]['mRNA'].split('\t')
                    if obj.seqids() == token[0]:
                        print_flag = False
                        ### print new mRNA
                        out.write(token[0]+ '\t'+ source + '\t' + token[2]+ '\t' + token[3]+ '\t'+ token[4]+'\t'+ token[5]+ '\t'+ token[6]+'\t'+ token[7]+ '\t'+ (token[8].split("Parent=")[0]).replace(reaplcement_IDs[ID], ID)+"Parent="+g_id+";Name="+ID+'\n')
                        
                        ### print exon lines of the lines
                        for i in hash_combine[reaplcement_IDs[ID]]['exon'].split(',')[1:]:
                            i = i.replace(reaplcement_IDs[ID], ID)
                            token = i.split('\t')
                            out.write(token[0]+'\t'+ \
                                      source + '\t' + \
                                      '\t'.join(token[2:])+'\n')
                            
                        ### print CDS lines of the lines
                        CDS_count = 0
                        for i in hash_combine[reaplcement_IDs[ID]]['CDS'].split(',')[1:]:
                            CDS_count += 1
                            i = i.replace(reaplcement_IDs[ID], ID)
                            token = i.split('\t')
                            out.write(token[0]+ '\t'+ source + '\t' + token[2]+ '\t' + token[3]+ '\t'+ token[4]+'\t'+ token[5]+ '\t'+ token[6]+'\t'+ token[7]+ '\t'+ (token[8].replace(reaplcement_IDs[ID], ID)).split("ID=")[0]+"ID="+ID+'.CDS.'+str(CDS_count)+";Parent="+ID+'\n')                    
                else:
                    print_flag = True
                    out.write(line+'\n')
            else:
                if print_flag == True:
                    out.write(line+'\n')
    out.close()
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### hash the replacement ids
    reaplcement_IDs, hash_combine_ID = hash_replacementIds()
    
    ### hash combine gff3 file with IDs to be replaced
    hash_combine = hash_combined(hash_combine_ID)
    
    ### print new GFF3 file
    make_gff3(reaplcement_IDs, hash_combine)
    
    ### close the logfile
    o.close()