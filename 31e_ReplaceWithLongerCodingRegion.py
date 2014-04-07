#-----------------------------------------------------------+
#                                                           |
# 31e_ReplaceWithLongerCodingRegion.py - Script to find the |
# longest protein coding evidence with overlapping exons    |
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
import os,sys,getopt, re, threading
from multiprocessing import Process, Queue
import classGene

### global variables
global last_gff3, evidences_gff3, gene_size_difference, exon_count_difference, min_exon_overlap, max_exonic_differences, min_coding_differences

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
            python 31e_ReplaceWithLongerCodingRegion.py
                                    -l/--last_gff3
                                    -e/--evidences_gff3
                                    -g/--gene_size_difference
                                    -c/--exon_count_difference
                                    -o/--min_exon_overlap
                                    -d/--max_exonic_differences
                                    -m/--min_coding_differences [in fraction, default 0.2]
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    
    global outfile
    #outfile = open(str(now.strftime("%Y-%m-%d_%H%M."))+'LongerCDS.txt','w')
    outfile = open(str("temp.")+'LongerCDS.txt','w')
    
    global last_gff3, evidences_gff3, gene_size_difference, exon_count_difference, min_exon_overlap, max_exonic_differences, min_coding_differences
    gene_size_difference=20
    min_coding_differences=0.2
    max_exonic_differences=50
    
    
    try:
        opts, args = getopt.getopt(argv,"hl:e:g:c:o:d:m:",\
                                   ["last_gff3=",\
                                    "evidences_gff3=",\
                                    "gene_size_difference=", \
                                    "exon_count_difference=",\
                                    "min_exon_overlap=",\
                                    "max_exonic_differences=",\
                                    "min_coding_differences="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-l", "--last_gff3"):
            last_gff3 = arg
        elif opt in ("-e", "--evidences_gff3"):
            evidences_gff3 = arg
        elif opt in ("-g", "--gene_size_difference"):
            gene_size_difference = arg
        elif opt in ("-c", "--exon_count_difference"):
            exon_count_difference = arg
        elif opt in ("-o", "--min_exon_overlap"):
            min_exon_overlap = arg
        elif opt in ("-d", "--max_exonic_differences"):
            max_exonic_differences = arg
        elif opt in ("-m", "--min_coding_differences"):
            min_coding_differences = arg
            

    
def get_size(file):
    count = 0
    chroHash = {}
    for line in open(file,'r'):
        count += 1
        ### print the lines processed
        if count%100000 == 0:
            print 'Lines processed: ', '{:9,.0f}'.format(count)
        line = line.strip()
        if len(line)>1 and not line.startswith('#'):
            token = line.split('\t')
            if token[0] not in chroHash:
                chroHash[token[0]] = ''
    return chroHash

                

def hash_last_gff3_CDS(chromosome):
    CDS = {}
    count = 0
    exons = {}
    for line in open(last_gff3,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            ### print the lines processed
            count += 1
            if count%100000 == 0:
                print 'Hashing CDS: ', '{:9,.0f}'.format(count)
            ### check if it is gff3 format line
            if re.search('ID=',line):
                obj = classGene.GFF3(line)
                if obj.types()=='mRNA' and obj.seqids() == chromosome:
                    ID = str(obj)
                    CDS[ID] = 0
                    exons[ID] = []
                if obj.types()=='exon' and obj.seqids() == chromosome:
                    exons[ID].append( (int(obj.starts()),int(obj.ends())) )
                if obj.types()=='CDS' and obj.seqids() == chromosome:
                    if ID in CDS:
                        CDS[ID] += (int(obj.ends()) - int(obj.starts()) + 1)                  
    
    return exons, CDS

def find_gene_overlaps(evidences, exons_transcript, CDS_transcript, exons, CDS, chromosome):
    count = 0
    for line in open(last_gff3,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            count += 1
            if count%10000 == 0:
                print 'find_gene_overlaps: ', '{:9,.0f}'.format(count)
            obj = classGene.GFF3(line)
            if obj.types()=='mRNA' and obj.seqids()==chromosome:
                ID = str(obj)
                ### find the best replacement
                delta_CDS = []
                evi_list = []
                larger_found = False
                for i in evidences:
                    if int(obj.starts()) in evidences[i] and int(obj.ends()) in evidences[i]:
                        if ID !=  evidences[i][int(obj.starts())] and ID != evidences[i][int(obj.ends())]:
                            evidence_id = evidences[i][int(obj.starts())]
                            if CDS[ID] + float(min_coding_differences)*CDS[ID] < CDS_transcript[evidence_id] and CDS[ID] >= 0:
                                ### make sure that two IDs are isoformorms
                                if '.'.join(ID.split('.')[:-1]) != '.'.join(evidence_id.split('.')[:-1]):
                                    ### check the exonic difference
                                    set_exons_transcript = {}
                                    for (start,end) in exons_transcript[evidence_id]:
                                        for j in range(start, end+1):
                                            set_exons_transcript[j] = ''
                                    set_exons_transcript=set(set_exons_transcript)
                                    set_exons = {}
                                    for start,end in exons[ID]:
                                        for j in range(start, end+1):
                                            set_exons[j] = ''
                                    set_exons=set(set_exons)
                                    if abs(int(len(set_exons_transcript - set_exons))) < int(max_exonic_differences):
                                        delta_CDS.append(CDS_transcript[evidence_id] - CDS[ID])
                                        evi_list.append(ID+'\t'+evidence_id+'\t'+i+'\t'+str(obj.seqids())+'\t'+str(obj.starts())+ '\t'+ str(CDS_transcript[evidence_id] - CDS[ID]) + '\t' + str(len(set_exons_transcript - set_exons)))
                                        larger_found = True
                if larger_found == True:
                    outfile.write(evi_list[delta_CDS.index(min(delta_CDS))]+'\n')
                


def hash_evidences(chromosome, exons, CDS):
    evidences = {} ### hash exonic co-ordinates by evidence
    CDS_transcript = {}
    exons_transcript = {}
    count = 0
    for line in open(evidences_gff3,'r'):
        line = line.strip()
        count += 1
        if len(line) > 0 and not line.startswith('#'):
            if count%100000 == 0:
                print 'Hashing evidence: ', '{:9,.0f}'.format(count)
            ### check if it is gff3 format line
            if re.search('ID=',line):
                obj = classGene.GFF3(line)
                if obj.types()=='mRNA' and obj.seqids() == chromosome:
                    ID = str(obj)
                    CDS_transcript[ID] = 0
                    exons_transcript[ID] = []
                    ### print the lines processed
                    if obj.sources() not in evidences:
                        evidences[obj.sources()] = {}
                    for i in range(int(obj.starts())-int(gene_size_difference)-1 ,int(obj.starts())+int(gene_size_difference)+1):
                        evidences[obj.sources()][i] = str(obj)
                    for i in range(int(obj.ends())-int(gene_size_difference)-1,int(obj.ends())+int(gene_size_difference)+1):
                        evidences[obj.sources()][i] = str(obj)
                if obj.types()=='exon' and obj.seqids() == chromosome:
                    exons_transcript[ID].append( (int(obj.starts()),int(obj.ends())) )
                if obj.types()=='CDS' and obj.seqids() == chromosome:
                    if ID in CDS_transcript:
                        CDS_transcript[ID] += (int(obj.ends()) - int(obj.starts()) + 1)
                   
    
    find_gene_overlaps(evidences, exons_transcript, CDS_transcript, exons, CDS, chromosome)

def parse(chromosome):
    
    exons, CDS = hash_last_gff3_CDS(chromosome)
    
    hash_evidences(chromosome, exons, CDS)



if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    print 'Hashing the chromosomes name'
    chroHash = get_size(last_gff3)
    
    for chromosome in sorted(chroHash):
        parse(chromosome)
    
    '''
    ### multithreading
    thread_list = []
    for chromosome in sorted(chroHash):
        t = Process(target=parse, args=(chromosome,))
        t.start()
        thread_list.append(t)
    for thread in thread_list:
        thread.join()
    ''' 
    
    ### close the logfile
    outfile.close()