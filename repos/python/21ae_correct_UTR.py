'''
Created on Feb 24, 2013
@author: vgupta
'''

# 21ae_correct_UTR.py - script to correct the UTR co-ordinates - /Users/vikas0633/Desktop/script/python
# Only parse sorted gff file to this script (order should be gene/mRNA/exon)

import os,sys,getopt, re
from C_loadFasta import *

### main argument to 

def options(argv):
    inputfile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print 'python 21ae_correct_UTR.py -i <inputfile>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python 21ae_correct_UTR.py -i <inputfile>'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
            
    return inputfile
    
                            
### a function to correct the UTRs co-ordinates
def correct_UTRs(elements,strand,chr,type,start,end,exon,CDS,mRNA_ID):
    if 'CDS' not in elements:
        print elements['exon'].strip()
    else:
        ### divide the correcting to strands
        print elements['exon'].strip()
        if strand == '+':
            ### hash the exon coordinates
            ### chr5    GeneMark.hmm    exon    24744063    24744193    .    +    .    ID=101603_t:exon:1;Name=101603_t:exon:1;Parent=101603_t
            ### get the lowest key in the exon
            min_key = sorted(exon)[0]
            flag = False
            flag2 = True
            for key in sorted(exon):
                start = exon[key][0]
                end = exon[key][1]
                if ((key == min_key) & (exon[key][1] < min(CDS)) & (exon[key][0] < min(CDS))): #print first exonic 5'UTR
                    flag = False
                    print chr +'\t'+ type + '\t' + 'five_prime_UTR_exon\t' + str(start) +'\t' + str(end) +'\t' + '.\t'+'+\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                
                elif ((exon[key][1] < min(CDS)) & (exon[key][0] < min(CDS))):
                    print chr +'\t'+ type + '\t' + 'five_prime_UTR_intron\t' + str(exon[key-1][1]+1) +'\t' + str(start-1) +'\t' + '.\t'+'+\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    print chr +'\t'+ type + '\t' + 'five_prime_UTR_exon\t' + str(start) +'\t' + str(end) +'\t' + '.\t'+'+\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'  
                
                elif ((exon[key][0] < min(CDS)) & (exon[key][1] > min(CDS)) & (flag == False)): #print last exonic 5'UTR
                    if key != min_key : ### do not print intron UTR if only one exon
                        print chr +'\t'+ type + '\t' + 'five_prime_UTR_intron\t' + str(exon[key-1][1]+1) +'\t' + str(start-1) +'\t' + '.\t'+'+\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    print chr +'\t'+ type + '\t' + 'five_prime_UTR_exon\t' + str(start) +'\t' + str(min(CDS)-1) +'\t' + '.\t'+'+\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    flag = True
                    flag2 = False
                
                if ((exon[key][1] > max(CDS)) & (flag2 == False)):
                    print chr +'\t'+ type + '\t' + 'three_prime_UTR_exon\t' + str(max(CDS)+1) +'\t' + str(end) +'\t' + '.\t'+'+\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    flag2 = True
                elif ((exon[key][1] > max(CDS)) & (flag2 == True)):
                    if key != min_key :
                        print chr +'\t'+ type + '\t' + 'three_prime_UTR_intron\t' + str(exon[key-1][1]+1) +'\t' + str(start-1) +'\t' + '.\t'+'+\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    print chr +'\t'+ type + '\t' + 'three_prime_UTR_exon\t' + str(start) +'\t' + str(end) +'\t' + '.\t'+'+\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
            print elements['CDS'].strip()
            
        if strand == '-':               ### process negative strand
            ### hash the exon coordinates
            ### chr5    GeneMark.hmm    exon    24744063    24744193    .    +    .    ID=101603_t:exon:1;Name=101603_t:exon:1;Parent=101603_t
            ### get the lowest key in the exon
            min_key = sorted(exon)[0]   ### exon with smallest co-ordinates 
            flag = False
            flag2 = True
            for key in sorted(exon):    ### go through each exon
                start = exon[key][0]
                end = exon[key][1]
                if ((key == min_key) & (exon[key][1] < min(CDS)) & (exon[key][0] < min(CDS))):          #print first exonic 5'UTR
                    flag = False
                    print chr +'\t'+ type + '\t' + 'three_prime_UTR_exon\t' + str(start) +'\t' + str(end) +'\t' + '.\t'+'-\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    
                elif ((exon[key][1] < min(CDS)) & (exon[key][0] < min(CDS))):
                    print chr +'\t'+ type + '\t' + 'three_prime_UTR_intron\t' + str(exon[key-1][1]+1) +'\t' + str(start-1) +'\t' + '.\t'+'-\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    print chr +'\t'+ type + '\t' + 'three_prime_UTR_exon\t' + str(start) +'\t' + str(end) +'\t' + '.\t'+'-\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'  
                    
                elif ((exon[key][0] < min(CDS)) & (exon[key][1] > min(CDS)) & (flag == False) & (exon[key][0] < min(CDS))): #print last exonic 5'UTR
                    if key != min_key : ### do not print intron UTR if only one exon
                        print chr +'\t'+ type + '\t' + 'three_prime_UTR_intron\t' + str(exon[key-1][1]+1) +'\t' + str(start-1) +'\t' + '.\t'+'-\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    print chr +'\t'+ type + '\t' + 'three_prime_UTR_exon\t' + str(start) +'\t' + str(min(CDS)-1) +'\t' + '.\t'+'-\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    flag = True
                    flag2 = False
                    
                if ((exon[key][1] > max(CDS)) & (flag2 == False)):
                    print chr +'\t'+ type + '\t' + 'five_prime_UTR_exon\t' + str(max(CDS)+1) +'\t' + str(end) +'\t' + '.\t'+'-\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    flag2 = True
                elif ((exon[key][1] > max(CDS)) & (flag2 == True)):
                    if key != min_key :
                        print chr +'\t'+ type + '\t' + 'five_prime_UTR_intron\t' + str(exon[key-1][1]+1) +'\t' + str(start-1) +'\t' + '.\t'+'-\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
                    print chr +'\t'+ type + '\t' + 'five_prime_UTR_exon\t' + str(start) +'\t' + str(end) +'\t' + '.\t'+'-\t.\t' +'ID='+mRNA_ID+'.'+str(key)+';'+'Parent='+mRNA_ID+';'
            print elements['CDS'].strip()            

def parse(file):
    elements = {}
    count = 0
    exon = {}
    CDS = []
    first_line = True
    for line in open(file,'r'):
        line = line.strip()
        token = line.split()
        if len(line) > 1:
            if line[0] != '#':
                line = line.replace('model.','gene.')
                line = line.replace('.mrna','.path')
                if (token[1] == "GeneMark.hmm") & (token[2] == "exon"):
                    line = line.replace("_t","_g")
                if token[1] == "TAU":
                    tokens=line.split('.')
                    line = '.'.join(tokens[0:len(tokens)-1])+';'    
                if token[2] == "gene":
                    if first_line == False:
                        correct_UTRs(elements,strand,chr,type,start,end,exon,CDS,mRNA) ### send the mRNA strand,start,end for UTR annotation
                        elements = {}
                        count = 0
                        exon = {}
                        CDS = []
                    first_line = False
                    print line
                if token[2] == "mRNA":
                    print line
                    if token[1]=="CUFFLINKS":
                        match = re.search(r'ID=.+;',line)
                        if match:
                            match = match.group().split(';')[0].replace('ID=','')
                            mRNA = match
                            match = match.replace('mrna','path')
                    else:
                        match = re.search(r'ID=.*',line)
                        if match:
                            match = match.group().split(';')[0].replace('ID=','')
                            mRNA = match
                            match = match.replace('mrna','path')
                    strand = token[6]
                    chr = token[0]
                    type = token[1]
                    start = int(token[3])
                    end = int(token[4])
                if (token[2] != "gene") & (token[2] != "mRNA"):
                    ### hash the exon/CDS co-ordinates
                    count += 1
                    if token[2]=="exon":
                        match_exon = re.search(r'Parent=.+',line)           ### make sure exons has correct mRNA parents
                        if match_exon:
                            match_exon = match_exon.group().split(';')[0].replace('Parent=','')
                            if match_exon != match:
                                line = line.split("Parent")[0]+'Parent='+match+';'
                        exon[count] = int(token[3]),int(token[4])
                    if token[2]=="CDS":
                        CDS.append(int(token[3]))
                        CDS.append(int(token[4]))
                    
                    
                    ### hash the Exon/CDS lines
                    if token[2] not in elements:
                        elements[token[2]] = line +'\n'
                    else:
                        elements[token[2]] += line +'\n'
    correct_UTRs(elements,strand,chr,type,start,end,exon,CDS,mRNA)

if __name__ == "__main__":
    
    file = options(sys.argv[1:])

    parse(file)
