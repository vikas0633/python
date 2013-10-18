#-----------------------------------------------------------+
#                                                           |
# 128_PFformatter.py - Script to create pathologic formatted file    |
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
global gff3, ec, go, fasta, HEADER

### Attributes
attributes = ['ID', 'NAME', 'STARTBASE', 'ENDBASE', 'PRODUCT-TYPE','SYNONYM','GENE-COMMENT','FUNCTION','EC','GO','DBLINK','//']

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
    global gff3, ec, go, fasta, HEADER
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hg:e:o:f:",["gff3=","EC=","GO=","fasta="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-g", "--gff3"):
            gff3 = arg
        elif opt in ("-e", "--EC"):
            ec = arg
        elif opt in ("-o", "--GO"):
            go = arg
        elif opt in ("-f", "--fasta"):
            fasta = arg
            
    logfile(gff3)
            
def hashEC():
    global gff3, ec, go, fasta,  HEADER
    '''
    s.id	q.id	q.length	s.length	identity	evalue	s.cov	protein	species	Pathway.id	Pathway.name	Reaction.id	EC	cyc
    834247-MONOMER	CUFF.177.1	4100	364	80.1	0.00E+00	99.5	caffeate 3-O-methyltransferase	Populus trichocarpa	PWY-2181	free phenylpropanoid acid biosynthesis	RXN-1104	EC-2.1.1.68	poplarcyc
    834247-MONOMER	CUFF.177.1	4100	364	80.1	0.00E+00	99.5	caffeate 3-O-methyltransferase	Populus trichocarpa	PWY-1121	suberin biosynthesis	RXN-1104	EC-2.1.1.68	poplarcyc
    AT1G02280-MONOMER	CUFF.27750.1	1411	298	62.1	6.00E-130	98.0	GTPase	Arabidopsis thaliana	NA	NA	NA	NA	NA
    '''
    t_id = '' ### transcript id - 2nd 
    EC = {} ### EC number - 13th
    Function = {} ### protein name - 8th
    Name = {} ### pathway name - 10th 
    GeneComment = {} ### ReactionID - 12th
    for line in open(ec, 'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            tokens = line.split('\t')
            t_id = tokens[1]
            if tokens[12] != 'NA':
                ec_no = tokens[12]
                if t_id in EC:
                    EC[t_id].append(ec_no)
                else:
                    EC[t_id] = [ec_no]
                    
                Function[t_id] = tokens[7]
                Name[t_id] = tokens[9]
                GeneComment[t_id] = tokens[11]
    
    return EC, Function, Name, GeneComment

def hashGO():
    '''
    Gene_ID	Transcript_ID	protein_ID	GO_ID	Category	Description
    10091_g	10091_g.1	10091_g.1_4_ORF1	GO:0003723	Molecular Function	RNA binding
    10091_g	10091_g.1	10091_g.1_4_ORF1	GO:0003964	Molecular Function	RNA-directed DNA polymerase activity
    10091_g	10091_g.1	10091_g.1_4_ORF1	GO:0006278	Biological Process	RNA-dependent DNA replication
    '''
    global gff3, ec, go,fasta,  HEADER
    
    t_id = '' ### 2nd
    ProductType = {} ### P
    GO = {}    ### 6th|4th | 22914220 | IMP 
    SYNONYM = {} ### 5th
    
    
    for line in open(go, 'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            tokens = line.split('\t')
            ProductType[t_id] = 'P'
            go_id = tokens[5]+'|'+ tokens[3][3:] + '|' + '22914220 | IMP'
            if t_id in GO:
                GO[t_id].append(go_id)
            else:
                GO[t_id] = [go_id]
                SYNONYM[t_id] = tokens[4]
    
    return ProductType, GO, SYNONYM

def parseGFF3(EC, Function, Name, GeneComment, ProductType, GO, SYNONYM, header, pf_out, g):
    ### Attributes
    #attributes = ['ID', 'NAME', 'STARTBASE', 'ENDBASE', 'PRODUCT-TYPE','SYNONYM','GENE-COMMENT','FUNCTION','EC','GO','DBLINK','//']

    for line in open(gff3, 'r'):
        line = line.strip()
        if len(line)>0 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            if obj.types() == 'mRNA' and obj.seqids() == header:
                t_id = str(obj)
                g.write('ID\t'+str(header)+'\n')
                g.write('Name\t'+str(header)+'\n')
                g.write('TYPE\t:CHRSM\n')
                g.write('CIRCULAR?\tN\n')
                g.write('ANNOT-FILE\t'+header+'.pf\n')
                g.write('SEQ-FILE\t'+header+'.fa\n')
                g.write('\\\\'+'\n')
                        

                if t_id in Name:
                    pf_out.write('ID\tVc_'+str(t_id)+'\n')
                    pf_out.write(attributes[1] +'\t'+ Name[t_id]+'\n')
                    pf_out.write(attributes[2] + '\t' + str(obj.starts())+'\n')
                    pf_out.write(attributes[3] + '\t' + str(obj.ends())+'\n')
                    if t_id in ProductType:
                        pf_out.write(attributes[4] + '\t' + str(ProductType[t_id])+'\n')
                    else:
                        pf_out.write(attributes[4] + '\t' + str('P')+'\n')
                    if t_id in SYNONYM:
                        pf_out.write(attributes[5] + '\t' + str(SYNONYM[t_id])+'\n')
                    if t_id in GeneComment:
                        pf_out.write(attributes[6] + '\t' + str(GeneComment[t_id])+'\n')
                    if t_id in Function:
                        pf_out.write(attributes[7] + '\t' + str(Function[t_id])+'\n')
                    if t_id in EC:
                        for j in range(len(EC[t_id])):
                            ec_no = EC[t_id][j].replace('EC-','')
                            if len(ec_no.split('.')) == 3:
                                ec_no +=  '.-'
                            pf_out.write(attributes[8] + '\t' + str(ec_no+'\n'))
                    if t_id in GO:
                        for j in range(len(GO[t_id])):
                            pf_out.write(attributes[9] + '\t' + str(GO[t_id][j])+'\n')
                            pf_out.write(attributes[10] + '\tGO:' + str(GO[t_id][j]).split('|')[1]+'\n')
                    pf_out.write(attributes[11]+'\n')

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    g = open('genetic-elements.dat','w')
    
    ### hash EC file
    EC, Function, Name, GeneComment = hashEC()
            
    ### hash GO file
    ProductType, GO, SYNONYM = hashGO()
    
    
    for line in open(fasta, 'r'):
        line = line.strip()
        if line.startswith('>'):
            header = line[1:].split()[0]
            print header
            fasta_out = open(header + '.fa', 'w')
            pf_out = open(header + '.pf', 'w')
            fasta_out.write('>'+header+'\n')
            
            ### print pathologic format file
            parseGFF3(EC, Function, Name, GeneComment, ProductType, GO, SYNONYM, header, pf_out, g)
        else:
            fasta_out.write(line+'\n')
    
    ### close the logfile
    o.close()