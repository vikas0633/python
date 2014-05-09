### 31c_TAU.py - script to add ORF to the gff3 file


'''
USAGES:
nice -n 19 python /Users/vikas0633/Desktop/script/python/21t_tau.py -f /Users/vikas0633/Desktop/temp/Ljchr0-6_pseudomol_20120830.chlo.mito.fa -g sample.20121108_merged_gene_models.gff3

'''
'''
custom script
It must do following:
1. Process one transcript at a time
2. Take out from targeted region from genome sequence
3. Replace the co-ordinates in the gff file according to the genome extract
4. Run TAU on extract genome/ transcript gene model
5. Add the CDS/UTRs on the gene structure
6. Remember +/- strand when adding ORF
7. Replace the co-ordinates back to original
8. Go back to step 1
'''

import os,sys, getopt
from D_longest_fasta_sequence_header import *

### GLOBAL VARIABLEs
global  infile, gff, intron_size
TAU_PATH="~/scripts/TAU"
SCRIPT_DIR="~/scripts"

#TAU_PATH="/u/vgupta/01_genome_annotation/tools/TAU"
#SCRIPT_DIR="/u/vgupta/script/python"

##TAU_PATH="~/Desktop/tools/TAU"
#SCRIPT_DIR="~/Desktop/script/python"


def file_empty(file):
    count = sum([1 for line in open(file)])
    if count == 0:
        sys.exit(file+' is empty')    



def options(argv):
    global  infile, gff, intron_size
    try:
        opts, args = getopt.getopt(argv,"hf:g:I:",["genome=","gff=","--intron_size="])
    except getopt.GetoptError:
        print 'python 21t_tau.py -f <genome> -g <gff>'
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print 'python 21t_tau.py -f <genome> -g <gff>'
            sys.exit()
        elif opt in ("-f", "--genome"):
            infile = arg
        elif opt in ("-g", "--gff"):
            gff = arg
        elif opt in ("-I", "--intron_size"):
            intron_size = arg 
    
def extract_seq(infile, header, mRNA_st, mRNA_en):
    os.system('echo ">"'+str(header)+' >temp.fa')
    os.system('nice -n 19 fastacmd -d '+infile+' -p F -s '+header+' -L '+str(mRNA_st)+','+str(mRNA_en) +'| awk "NR>1" >>temp.fa')

def exon(count,mRNA_id,mRNA_st,mRNA_en,strand,model_no,chr_id,no):
    exons = {}
    all_exons = {}
    exon_no = 0
    all_exon_no = 0
    file = './temp'+str(count)+'/temp'+str(count)+'.gff'
    for line in open(file,'r'):
        token = line.strip().split('\t')
        token=line.split('\t')
        start = int(token[3])
        end = int(token[4])
        if int(line.split(';')[0].split('.')[-2]) == no: ### check if same frame as longest cdna
            if token[2].startswith('CDS') & (chr_id == line.split(';')[0].split(' ')[-1]) :
                exon_no += 1
                #Ljchr1_pseudomol_20120830    TAU    CDSf    0    861    .    +    .    temp11 Ljchr1_pseudomol_20120830.1.1;evidence ORF;confidence 3.11627906976744;average_depth 1;
                ID = 'ID='+mRNA_id+'.'+str(exon_no)+';'+'Parent='+mRNA_id+'.exon.'+str(no); 
                #lin = token[0]+'\t'+token[1]+'\t'+'exon'+'\t'+str(start+mRNA_st)+'\t'+str(end+mRNA_st)+'\t'+token[5]+'\t'+'.'+'\t'+token[7]+'\t'+ID
                #o.write(lin+'\n')
                exons[exon_no] = start, end
        if token[2][0:3]=="CDS":
            all_exon_no += 1
            all_exons[all_exon_no] = start, end
    
    return exons, all_exons

def call_CDS(count,mRNA_id,gene_model,mRNA_st,mRNA_en,strand,model_no, exons, chr_id,no, all_exons):
    file = './temp'+str(count)+'/temp'+str(count)+'.cds.fa'
    CDS = False
    mRNA_strand = strand
    cds_strand = strand
    for line in open(file,'r'):
        if line.startswith('>'):
            if int(line.split(':')[0].split('.')[-2]) == no:
                # Ljchr1_pseudomol_20120830    Augustus    mRNA    21237    22777    .    +    .    ID=model.g28684.t1;Parent=gene.g28684.t1
                # >Ljchr1_pseudomol_20120830.1.1:cds:+:2366:3503
                token = line.split('.')
                tokens = line.split(':')
                start = int(tokens[3])
                end = int(tokens[4])
                cds_strand = tokens[2]
		
                ## changing CDS coordinating 2014-03-11
                mRNA_st += -1
                mRNA_en += -1
                #print cds_strand, strand, tokens[0][1:], chr_id
                if (chr_id == tokens[0][1:]):
                    print 'all_exons', all_exons
                    ID = 'ID='+mRNA_id+'.'+str(model_no)+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                    frame = 0
                    if cds_strand == '+':
                        CDS = True
                        #mRNA_st = mRNA_st + exons[1][0]
                        if len(exons)==1:
                            #five_prime_UTR_st =  exons[1][0] + mRNA_st -1
                            five_prime_UTR_st = mRNA_st -1
                            five_prime_UTR_en = start - exons[1][0] + mRNA_st -2
                            ID = 'ID='+mRNA_id+'.'+str('5_prime_UTR')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            lin = token[0][1:]+'\t'+'TAU'+'\t'+'five_prime_UTR'+'\t'+str(five_prime_UTR_st)+'\t'+str(five_prime_UTR_en)+'\t'+'.'+'\t'+cds_strand+'\t'+'.'+'\t'+ID
                            o.write(lin+'\n')
                            CDS_st = start
                            CDS_en = end
                            ID = 'ID='+mRNA_id+'.'+str('CDS.1')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(CDS_st + mRNA_st )+'\t'+str(CDS_en + mRNA_st -1 )+'\t'+'.'+'\t'+cds_strand+'\t'+str(frame)+'\t'+ID
                            o.write(lin+'\n')
                            cds_len = CDS_en - CDS_st + 1
                            exon_st = exons[1][0]
                            exon_en = exons[1][1]
                            last_cds_len = 0
                            three_prime_UTR_st = CDS_en + mRNA_st
                            
                            three_prime_UTR_en = exons[len(exons)][1] + mRNA_st - exons[1][0]
                            ID = 'ID='+mRNA_id+'.'+str('3_prime_UTR')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            if (three_prime_UTR_en - three_prime_UTR_st > 1):
                                lin = token[0][1:]+'\t'+'TAU'+'\t'+'three_prime_UTR'+'\t'+str(three_prime_UTR_st)+'\t'+str(three_prime_UTR_en)+'\t'+'.'+'\t'+cds_strand+'\t'+'.'+'\t'+ID
                                o.write(lin+'\n')
                                
                        elif (cds_strand == mRNA_strand):
                            exon_length = 0
                            introns = 0
                            for i in range(1,len(exons)+1):
                                exon_length += exons[i][1] - exons[i][0] +1
                                if i > 1:
                                    introns += exons[i][0] - last_exon -1
                                if exon_length >= start:
                                    j = i
                                    break
                                last_exon =  exons[i][1]
        
                            
                            five_prime_UTR_st = mRNA_st - exons[1][0]
                            five_prime_UTR_en = start + introns + mRNA_st -1
                            ID = 'ID='+mRNA_id+'.'+str('5_prime_UTR')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            lin = token[0][1:]+'\t'+'TAU'+'\t'+'five_prime_UTR'+'\t'+str(five_prime_UTR_st)+'\t'+str(five_prime_UTR_en)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
                            o.write(lin+'\n')
                            CDS_st = start + introns -1
                            if exons[j][1] - exons[1][0] - introns >= end:
                                if (start < (exons[1][1] - exons[1][0]))&(end > (exons[1][1] - exons[1][0])):
                                    lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str((start)+mRNA_st)+'\t'+str(exons[1][1]+ mRNA_st )+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
                                    o.write(lin+'\n')
                                three_prime_UTR_st = (end+introns) + 1 + mRNA_st
                                CDS_en = (end+introns) + 1 
                                flag_end = True
                            else:
                                CDS_en = exons[j][1] - exons[1][0] -1
                                flag_end = False
                            ID = 'ID='+mRNA_id+'.'+str('CDS.1')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(CDS_st+ mRNA_st + 1)+'\t'+str(CDS_en+ mRNA_st )+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
                            o.write(lin+'\n')
                            three_prime_UTR_st = CDS_en + mRNA_st + 1 - exons[1][0]
                            cds_len = CDS_en - CDS_st + 1
                            frame = 3 - (cds_len) % 3
                            if frame == 3:
                                frame = 0
                            if flag_end == False:
                                for i in range(j+1,len(exons)+1):
                                    exon_st = exons[i][0] - exons[1][0] -1
                                    exon_en = exons[i][1] - exons[1][0] -1
                                    last_cds_len = cds_len
                                    cds_len += (exon_en - exon_st) + 1
                                    ID = 'ID='+mRNA_id+'.'+str('CDS.')+str(i)+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            
                                    if (start <= exon_st) & (end-start > cds_len):
                                        lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(exon_st+ mRNA_st)+'\t'+str(exon_en+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
                                        o.write(lin+'\n')
                                        three_prime_UTR_st = exon_st+ mRNA_st + 1
                                        frame = 3 - (cds_len)%3
                                        if frame == 3:
                                            frame = 0
                                    elif end-start <= cds_len:
                                        lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(exon_st+ mRNA_st)+'\t'+str(exon_st+(end - start - last_cds_len)+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
                                        o.write(lin+'\n')
                                        three_prime_UTR_st = exon_st+(end - start - last_cds_len) + 3 + mRNA_st
                                        break
        
                            three_prime_UTR_en = exons[len(exons)][1] + mRNA_st - exons[1][0]
                            ID = 'ID='+mRNA_id+'.'+str('3_prime_UTR')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            if (three_prime_UTR_en - three_prime_UTR_st > 1):
                                lin = token[0][1:]+'\t'+'TAU'+'\t'+'three_prime_UTR'+'\t'+str(three_prime_UTR_st)+'\t'+str(three_prime_UTR_en)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
                                o.write(lin+'\n')
                    
                    if cds_strand == '-':
                        #mRNA_st = mRNA_st - exons[1][0]
                        CDS = True
                        if len(exons)==1:
                            five_prime_UTR_st =  mRNA_en
                            five_prime_UTR_en =  exons[len(exons)][1] - start  + mRNA_st - 1
                            ID = 'ID='+mRNA_id+'.'+str('5_prime_UTR')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            lin = token[0][1:]+'\t'+'TAU'+'\t'+'five_prime_UTR'+'\t'+str(five_prime_UTR_en)+'\t'+str(five_prime_UTR_st)+'\t'+'.'+'\t'+cds_strand+'\t'+'.'+'\t'+ID
                            o.write(lin+'\n')
                            ID = 'ID='+mRNA_id+'.'+str('CDS.1')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            if len(all_exons) == 1:
                                CDS_st = exons[len(exons)][1] - start -2 - exons[1][0] 
                                CDS_en = exons[len(exons)][1] - end -1 - exons[1][0]
                                
                            elif exons[1][1] <= all_exons[1][1]: ## single cds but in first exon
                                CDS_st = exons[len(exons)][1] - start -2 - exons[1][0] 
                                CDS_en = exons[len(exons)][1] - end -1 - exons[1][0]
                            
                            else:
                                CDS_st = exons[len(exons)][1] - start -1 - exons[1][0] 
                                CDS_en = exons[len(exons)][1] - end - exons[1][0]
                            
                            lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(CDS_en + mRNA_st +1 )+'\t'+str(CDS_st + mRNA_st + 1)+'\t'+'.'+'\t'+cds_strand+'\t'+str(frame)+'\t'+ID
                            o.write(lin+'\n')
                            cds_len = CDS_st - CDS_en + 1
                            frame = 3- (cds_len) % 3
                            if frame == 3:
                                frame = 0
                            last_cds_len = 0
                            three_prime_UTR_st = CDS_en - 1 + mRNA_st
                            three_prime_UTR_en = mRNA_st - 1
                            ID = 'ID='+mRNA_id+'.'+str('3_prime_UTR')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            lin = token[0][1:]+'\t'+'TAU'+'\t'+'three_prime_UTR'+'\t'+str(three_prime_UTR_en)+'\t'+str(three_prime_UTR_st)+'\t'+'.'+'\t'+cds_strand+'\t'+'.'+'\t'+ID
                            o.write(lin+'\n')
                            
                        elif (cds_strand == mRNA_strand):
                            exon_length = 0
                            introns = 0
                            for i in range(len(exons),0,-1):
                                exon_length += exons[i][1]-exons[i][0] +1
                                if i < len(exons):
                                    introns += last_exon - exons[i][1] -1
                                if exon_length > start:
                                    j = i
                                    break
                                last_exon =  exons[i][0]
                            five_prime_UTR_st =  mRNA_en
                            five_prime_UTR_en =  exons[len(exons)][1] - start - introns + mRNA_st - 1
                            ID = 'ID='+mRNA_id+'.'+str('5_prime_UTR')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            lin = token[0][1:]+'\t'+'TAU'+'\t'+'five_prime_UTR'+'\t'+str(five_prime_UTR_en)+'\t'+str(five_prime_UTR_st)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
                            o.write(lin+'\n')
                            ID = 'ID='+mRNA_id+'.'+str('CDS.1')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            
                               
                            CDS_st = exons[len(exons)][1] - start - introns -2
                            
                            
                            if exons[len(exons)][1] - exons[j][0] - introns >= end:
                                CDS_en = exons[len(exons)][1] - (end+introns) -2
                                cds_len = CDS_st - CDS_en + 1
                                frame =  (cds_len) % 3 -1
                                if frame == 3:
                                    frame=0

                                if (start < (exons[len(exons)][1] - exons[len(exons)][0])) & (end > (exons[len(exons)][1] - exons[len(exons)][0])):
                                    lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(exons[len(exons)][0]+ mRNA_st +1 )+'\t'+str(exons[j][1] - start + mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
                                    o.write(lin+'\n')
                                three_prime_UTR_st = exons[len(exons)][1] - (end+introns) - 1 + mRNA_st
                                flag_end = True
                            else:
                                CDS_en = exons[j][0] -1
                                flag_end = False
                            cds_len = CDS_st - CDS_en + 1
                            frame =  (cds_len) % 3
                            if frame == 3:
                                frame = 0
                            
                            print CDS_en, CDS_st
                            if CDS_en <= CDS_st + 3 :
                                lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(CDS_en+ mRNA_st)+'\t'+str(CDS_st+ mRNA_st + 1)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
                                o.write(lin+'\n')
                            three_prime_UTR_st = CDS_en + mRNA_st -1
                            if flag_end == False:
                                for i in range(j-1,0,-1):
                                    exon_st = exons[i][0] -1 
                                    exon_en = exons[i][1] -1
                                    last_cds_len = cds_len
                                    cds_len += (exon_en - exon_st) + 1
                                    ID = 'ID='+mRNA_id+'.'+str('CDS.')+str(i)+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                                    #print start, exon_st, end-start, cds_len
                                    if (start <= exon_st) & (end-start > cds_len):
                                        frame = (cds_len) % 3
                                        if frame == 3:
                                            frame = 0
                                        lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(exon_st+ mRNA_st)+'\t'+str(exon_en+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
                                        o.write(lin+'\n')
                                        three_prime_UTR_st = exon_st+ mRNA_st -1
                                    elif end-start <= cds_len:
                                        lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(exons[i][1]-(end -start -last_cds_len)+mRNA_st + 1)+'\t'+str(exons[i][1]+ mRNA_st -1) +'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
                                        o.write(lin+'\n')
                                        three_prime_UTR_st = exons[i][1]-(end -start -last_cds_len) - 1 + mRNA_st
                                        break
                            three_prime_UTR_en = mRNA_st - 1
                            ID = 'ID='+mRNA_id+'.'+str('3_prime_UTR')+';'+'Parent='+mRNA_id+'.'+str(no)+';'
                            lin = token[0][1:]+'\t'+'TAU'+'\t'+'three_prime_UTR'+'\t'+str(three_prime_UTR_en)+'\t'+str(three_prime_UTR_st)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
                            o.write(lin+'\n')
    return cds_strand
                
            
def mRNA(count,gene_id,gene_model,mRNA_st,mRNA_en,strand,model_no):
    file = './temp'+str(count)+'/temp'+str(count)+'.cdna.fa'
    file_cds = './temp'+str(count)+'/temp'+str(count)+'.cds.fa'
    Longheader = longest_seq(file_cds)
    flag = False ### flag for print mRNA co-ordinates
    for line in open(file,'r'):
        line = line.strip()
        if line.startswith('>'):
            no = int(line.split(':')[0].split('.')[-2])
            token = line.split('.')
            tokens = line.split(':')
            start = int(tokens[2]) + mRNA_st
            end = int(tokens[3]) + mRNA_st
            ID = 'ID='+gene_id+'.'+str(model_no)+';'+'Parent='+gene_id+';'
            mRNA_id = gene_id+'.'+str(model_no)
            #strand = tokens[4].strip()                                                 # Do not take strand from cdna file
            #lin = token[0][1:]+'\t'+gene_model+'\t'+'mRNA'+'\t'+str(start)+'\t'+str(end)+'\t'+'.'+'\t'+tokens[4].strip()+'\t'+'.'+'\t'+ID
            lin = token[0][1:]+'\t'+gene_model+'\t'+'mRNA'+'\t'+str(mRNA_st)+'\t'+str(mRNA_en)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
            cds_strand = strand
            if flag == False:
                o.write(lin+'\n')
                flag = True
            if (Longheader.split(':')[0] == line.split(':')[0]) and (Longheader.split(':')[2]==line.split(':')[4]):            
                # Ljchr1_pseudomol_20120830    Augustus    mRNA    21237    22777    .    +    .    ID=model.g28684.t1;Parent=gene.g28684.t1
                # Ljchr1_pseudomol_20120830.1.1:cdna:1:1541:+
                mRNA_st = start -1
                ### check if the longest cdna is also in same strand
                #if strand == tokens[4].strip():
                ### call exons from TAU gff file
                exons, all_exons = exon(count,mRNA_id,start,end,strand,model_no,tokens[0][1:],no)
                ### call CDS/UTRs
                cds_strand = call_CDS(count,mRNA_id,gene_model,start,end,strand,model_no,exons,tokens[0][1:],no, all_exons)
                
                print Longheader, cds_strand
                
    return Longheader.split(':')[0], cds_strand

def process_transcript(count,gene_id,gene_model,mRNA_st,mRNA_en,strand,model_no):
    ### open cdna file and o.write( as mRNA gene-model
    return mRNA(count,gene_id,gene_model,mRNA_st,mRNA_en,strand,model_no)
    
def write_fasta(ID, infile, outfile, header, cds_strand):
    printed = False
    flag = False
    o = open(outfile,'a')
    for line in open(infile,'r'):
        if line.startswith('>'):
            if flag == True:
                    break
            if line.startswith(header):
                if cds_strand in line:
                    o.write('>'+ID+'\n')
                    printed = True
                    flag = True
               
        elif flag==True:
            o.write(line)
            
            
    if printed == False:
        for line in open(infile,'r'):
            if line.startswith('>'):
                o.write('>'+ID+'\n')
                if flag == True:
                    break
                else:
                    flag = True
            elif flag == True:
                o.write(line)
    o.close()
    

def save_files(ID, count, header, cds_strand):
    ### save the files for cdna, cds and proteins
    ### make a function to write the output file
    infile = './temp'+str(count)+'/temp'+str(count)+'.cdna.fa'
    outfile = 'TAU_cdna.fa'
    write_fasta(ID, infile, outfile, header, cds_strand)
    
    infile = './temp'+str(count)+'/temp'+str(count)+'.cds.fa'
    outfile = 'TAU_cds.fa'
    write_fasta(ID, infile, outfile, header, cds_strand)
    
    infile = './temp'+str(count)+'/temp'+str(count)+'.protein.fa'
    outfile = 'TAU_protein.fa'
    write_fasta(ID, infile, outfile, header, cds_strand)
    
    infile = './temp'+str(count)+'/temp'+str(count)+'.introns.fa'
    outfile = 'TAU_introns.fa'
    write_fasta(ID, infile, outfile, header, cds_strand)
    

def process_gff(gff,infile):
    id = ''
    count = -1
    gene_flag = ''
    new_transcript = False
    for line in open(gff,'r'):
        if len(line) >1:
            if line[0] != '#':
                line = line.strip()
                token = line.split('\t')
                if token[2] == 'gene':
                    gene_name = line
                    gene_model = token[1]
                    gene_id = (line.split('=')[1]).split(';')[0]
                    if gene_id != gene_flag:
                        gene_flag = gene_id
                        model_no = 1
                        gene_name = gene_name.replace('+','.').replace('-','.')
                        o.write(gene_name+'\n')
                    else:
                        gene_name = gene_name.replace('+','.').replace('-','.')
                        o.write(gene_name+'\n')
                if (token[2] == 'mRNA'):
                    count += 1 
                    # Ljchr1_pseudomol_20120830    Augustus    mRNA    21237    22777    .    +    .    ID=model.g28684.t1;Parent=gene.g28684.t1
                    if new_transcript == True:
                        tau_in.close()
                        ### make genome fasta file
                        extract_seq(infile, header, mRNA_st, mRNA_en)
                        ### change the format of the gff3 file
                        os.system('nice -n 19 python '+SCRIPT_DIR+'/21u_make_gff2.py -i tau_in -o temp.gff')
                        ### run TAU
                        os.system('nice -n 19 perl '+TAU_PATH +'/TAU.pl -I '+ str(intron_size)+' -A temp.fa -G temp.gff -O temp'+str(count) +' ')    
                        file = './temp'+str(count)+'/temp'+str(count)+'.cdna.fa'
                        
                        
                        if os.path.isfile(file):
                            if sum([1 for line1 in open(file)]) != 0:
                                ### process transcript
                                hedaer, cds_strand = process_transcript(count,last_gene_id,last_gene_model,mRNA_st,mRNA_en,strand,last_model_no)
                                save_files(last_gene_id+'.'+str(last_model_no),count,hedaer, cds_strand)
                                o.write(exon_line)                
                        
                        #sys.exit(0)
                    new_transcript = True
                    header = token[0]
                    mRNA_st = int(token[3]) 
                    mRNA_en = int(token[4]) 
                    strand = token[6]
                    start = int(token[3]) - mRNA_st +1 ### adding plus 1 on 2014-03-11
                    end = int(token[4]) - mRNA_st
                    lin = token[0]+'\t'+token[1]+'\t'+token[2]+'\t'+str(start)+'\t'+str(end)+'\t'+token[5]+'\t'+token[6]+'\t'+token[7]+'\t'+token[8]
                    tau_in = open('tau_in','w')
                    id = line.split('ID=')[1].split(';')[0]
                    tau_in.write(lin+'\n')
                    if model_no == 1:
                        last_gene_id = gene_id
                        last_gene_model = gene_model
                    last_model_no = model_no
                    model_no += 1
                    exon_line = ''
                    
                if (token[2] == 'exon'):
                    if line.split('Parent=')[1].split(';')[0] == id:
                        start = int(token[3]) - mRNA_st + 1 ### adding plus 1 on 2014-03-11
                        end = int(token[4]) - mRNA_st + 1
                        lin = token[0]+'\t'+token[1]+'\t'+token[2]+'\t'+str(start)+'\t'+str(end)+'\t'+token[5]+'\t'+token[6]+'\t'+token[7]+'\t'+token[8]
                        tau_in.write(lin+'\n')
                        exon_line += line+'\n'
        
        
    tau_in.close()
    count += 1
    model_no += 1
    if model_no == 1:
        gene_name = gene_name.replace('+','.').replace('-','.')
        o.write(gene_name+'\n')
    ### make genome fasta file
    extract_seq(infile, header, mRNA_st, mRNA_en)
    ### change the format of the gff3 file
    os.system('nice -n 19 python '+SCRIPT_DIR+'/21u_make_gff2.py -i tau_in -o temp.gff')

    ### run TAU
    os.system('nice -n 19 perl '+TAU_PATH +'/TAU.pl -I '+ str(intron_size)+' -A temp.fa -G temp.gff -O temp'+str(count) +' ')  
    
    
    file = './temp'+str(count)+'/temp'+str(count)+'.cdna.fa'
    if sum([1 for line in open(file)]) != 0:
        ### process transcript
        hedaer, cds_strand = process_transcript(count,last_gene_id,last_gene_model,mRNA_st,mRNA_en,strand,last_model_no)
        save_files(last_gene_id+'.'+str(last_model_no), count, hedaer, cds_strand)
        o.write(exon_line)
                
                
def remove_temp_files():
    os.system('[ -f tau_in ] && rm -r TAU_* || echo "No trash in the directory"')
    os.system('[ -f tau_in ] && rm -r temp* || echo "No trash in the directory"')

if __name__ == "__main__":

    options(sys.argv[1:])
    
    ### sort the gff file according to gene name
    hash_gff_gene = {}
    hash_gff_start = {}
    hash_gff_mRNA = {}
    hash_gff_exon = {}
    mRNA_id = {}
    
    ### remove TAU output file if exits
    remove_temp_files()
    
    ### format ref database
    os.system('formatdb -i '+infile+' -p F -o T')
    
    for line in open(gff,'r'):
        line = line.strip()
        token = line.split('\t')
        if len(line)>1:
            if line[0] != '#':
                if token[2] == "gene":
                    if (token[1]== "CUFFLINKS") or (token[1]=="Cufflinks"):
                        id = line.split('Name=')[1].split('%')[-1][2:]
                        hash_gff_gene[id]=line
                    else:
                        id = line.split('ID=')[1].split(';')[0]
                        hash_gff_gene[id]=line
                if token[2] == "mRNA":
                        if (token[1]== "CUFFLINKS") or (token[1]=="Cufflinks"):
                            id = line.split('ID=')[1].split(';')[0]
                            hash_gff_mRNA[id] = line
                            mRNA_id[id] = id 
                        else:
                            id = line.split('ID=')[1].split(';')[0]
                            parent_id = line.split('Parent=')[1].split(';')[0]
                            hash_gff_mRNA[parent_id] = line
                            mRNA_id[parent_id] = id
                if token[2] == "exon":
                    parent_id = line.split('Parent=')[1].split(';')[0]
                    id = line.split('ID=')[1].split(';')[0]
                    if parent_id in hash_gff_exon:
                        hash_gff_exon[parent_id] += '\n'+line 
                    else:
                        hash_gff_exon[parent_id] =line
    
    o = open(gff+'.out','w')
    for key in sorted(hash_gff_gene.keys()):
        gene_id = key
        gene = hash_gff_gene[key]
        o.write(gene+'\n')
        o.write(hash_gff_mRNA[gene_id]+'\n')
        o.write(hash_gff_exon[mRNA_id[gene_id]]+'\n')
    o.close()
    print "finished shorting the gff file"
    gff = gff+'.out'
    o = open('TAU_genemodel.gff3','w')
    ### process gene model file
    process_gff(gff,infile)
    
    os.system('rm -r temp*')
