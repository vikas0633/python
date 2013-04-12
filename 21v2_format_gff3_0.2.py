#21v2_format_gff3.py- script to format gff3 file in order to put in MySQL table sequal to 21v - /Users/vikas0633/Desktop/script/python

# Usage: python /Users/vikas0633/Desktop/script/python/21v2_format_gff3.py -i sample.gff3
'''
Usage: python /u/vgupta/script/python/21v2_format_gff3.py \
-i 20121227_conserved_Ljr3_0.refined_mRNA.sorted.gff3 \
-g Ljr_cdna.chr0-6.20130102.refined.fa \
-c Ljr_cds.chr0-6.20130102.refined.fa \
-p Ljr_cds_protein.Ljr3.0.20130102.refined.fa \
-b search_summary_length.txt \
-r 20121227_conserved_Ljr3_0.refined.ExonicRepeats \
-o mclGroups.20130102.txt \
-f orthoMCLfamilies.txt \
-s Ljr_cds_protein.Ljr3.0.20130102.refined.fa.blastout.genewise \
-j /u/vgupta/01_genome_annotation/23_ljrep_blast/Ljr_cds.chr0-6.20130102.refined.fa.blastout \
> 20130103_Ljr3_0_genemodel

																
'''

import os,sys,getopt, re


### main argument to 

def options(argv):
    inputfile = ''
    genomefile = ''
    cdsfile = ''
    proteinfile = ''
    homologyfile= ''
    blastfile = ''
    oMCLfamily = ''
    ljrep_blast = ''
    try:
       opts, args = getopt.getopt(argv,"hi:g:c:p:b:r:o:f:s:j:",["ifile=","gfile=","cfile=","pfile=", "bfile=","rfile=","orthoMCL=","oMCLfamily=","blast_search=","ljr_repeat="])
    except getopt.GetoptError:
       print 'python 21v_format_gff3.py -i <inputfile> -g <gfile> -c <cfile> -p <pfile> -b  <homology> -r <rfile> -o <orthoMCL> -f <oMCLfamily> -s <blast_search> -j <ljr_repeat>' 
       sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print 'python 21v_format_gff3.py -i <inputfile> -g <gfile> -c <cfile> -p <pfile> -b  <homology> -r <rfile> -o <orthoMCL> -f <oMCLfamily> -s <blast_search> -j <ljr_repeat>'
    		sys.exit()
    	elif opt in ("-i", "--ifile"):
    		inputfile = arg
        elif opt in ("-g", "--gfile"):
        	genomefile = arg
        elif opt in ("-c", "--cfile"):
        	cdsfile = arg
        elif opt in ("-p", "--pfile"):
        	proteinfile = arg
        elif opt in ("-b", "--homology"):
        	homologyfile = arg
        elif opt in ("-r", "--rfile"):
        	exon_repeat = arg
        elif opt in ("-o", "--orthoMCL"):
        	oMCLfile = arg
        elif opt in ("-s", "--blast_search"):
        	blastfile = arg
        elif opt in ("-f", "--oMCLfamily"):
        	oMCLfamily = arg
        elif opt in ("-j", "--ljr_repeat"):
        	ljrep_blast = arg
        
    return inputfile, genomefile, cdsfile, proteinfile, homologyfile, exon_repeat, oMCLfile, oMCLfamily, blastfile, ljrep_blast

def hashFasta(file):
	first_line = True
	seqs = {}
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0 :			
			if line[0] == '>':
				if first_line == False:
					seqs[header] = seq
				seq = ''
				header = line[1:]
			else:
				seq += line
		first_line = False			
	seqs[header] = seq

	return seqs
		
def exonRepeat(exon_repeat):
	exons = {}
	header = True
	for line in open(exon_repeat,'r'):
		line = line.strip()
		token = line.split('\t')
		if header == False:
			exons[token[0]] = str(token[1])+'\t'+str(token[2])+'\t'+str(token[4])
		header = False
	return exons
	
def hashHomology(hfile):
	homology = {}
	for line in open(hfile,'r'):
		line = line.strip()
		token = line.split('\t')
		
		ID = token[0]
		
		Mt_gene = token[10]
		Soybean_gene = token[13]
		
		homology[ID] = Mt_gene+'\t'+Soybean_gene
		
	return homology


def hashOrthoMCL(oMCLfile):
	orthoMCL = {}
	for line in open(oMCLfile,'r'):
		line = line.strip()
		token = line.split('\t')
		orthoMCL[token[0]] = token[1]
	return orthoMCL
	
	
def hashOrthofamily(oMCLfamily):
	orthofamily = {}
	for line in open(oMCLfamily,'r'):
		line = line.strip()
		token = line.split('\t')
		orthofamily[token[0]] = token[1] +'\t' +token[2] +'\t' +token[3] +'\t' +token[4] +'\t' +token[5] +'\t' +token[6] +'\t' + token[7] 
	return orthofamily

def hashblast(blastfile):
	blast = {}
	for line in open(blastfile,'r'):
		line = line.strip()
		token = line.split('\t')
		blast[token[0]] = token[1] +'\t'+ token[2] +'\t'+ token[3]
	return blast

def hashljrblast(file):
	''' This function is for hashing the blast output of Lj repeat '''
	blast = {}
	for line in open(file,'r'):
		line = line.strip()
		token = line.split('\t')
		blast[token[0]] = token[1] +'\t'+ token[2] +'\t'+ token[10]
	return blast

def format(inf, cDNA, cds, protein, homology, exons, orthoMCL, orthofamily, blast, ljrblast):
	count = 0
	### print header
	print 'Number'+'\t'+'Source'+'\t'+'ID'+'\t'+'Chromosome'+'\t'+'Start'+'\t'+'End'+'\t'+'Strand'+'\t'+'cDNASeq'+'\t'+'cDNASeq_length'+'\t'+'cdsSeq'+'\t'+'cdsSeq_len'+'\t'+'proteinSeq'+'\t'+'proteinSeq_len'+'\t'+'BlastHit[Medicago]'+'\t'+'BlastHit[Soyabean]'+'\t'+'ExonCount'+'\t'+'ExonicLength'+'\t'+'RepeatFraction' + '\t' + 'CDS/mRNA length ratio' + '\t' + 'OrthoMCL group' + '\t' +'ljr\t'+'mtr\t'+'gmax\t'+'ath\t'+'ptr\t'+'osa\t'+'sob\t'+ 'Annotation [UniProtBlast]' +'\t'+ '%Identity_UniProt' +'\t'+'Evalue_UniProt' + '\t' + 'Protein [Complete/Partial]' + '\tgene_name' + '\t' + 'Annotation [LjRep]' +'\t'+ '%Identity_LjRep' +'\t'+'Evalue_LjRep'
	
	for line in open(inf,'r'):
		token = line.split('\t')
		if len(token) > 3:
			if token[2] == "mRNA":
				count += 1
				line = line.strip()
				token = line.split('\t')
				if token[1] == "Augustus":
					source =  "Augustus"
				elif token[1] == "CUFFLINKS":
					source =  "Cufflinks"
				elif token[1] == "GeneMark.hmm":
					source =  "GeneMark"
				elif token[1] == "GlimmerHMM":
					source =  "GlimmerHMM"
				elif token[1] == "Ljchr0-6_pseudomol_20120830.chlo.mito.fa":
					source =  "Velvet"
				elif token[1] == "Ljchr30":
					source =  "Experimentally validated [NCBI]"
				elif token[1] == "Lotus_TC":
					source =  "Lotus_TC"
					
				### get the mRNA ID
				match = re.search(r'ID=.+;',line)
				if match:
					match = match.group().split(';')[0].replace('ID=','')
					ID = match				
				
				### get the gene name
				### idea is to split the mRNA name using '.' and the concatenate all back except the last one
				dot = ID.split('.')
				gene_name = dot[0]
				for i in range(1,len(dot)-1):
					gene_name += '.'+dot[i]
				
				
				### use the cdna seq
				cdnaSeq = ' '
				if ID in cDNA:
					cdnaSeq = cDNA[ID]
				
				cdsSeq = ''
				if ID in cds:
					cdsSeq = cds[ID]
				
				proteinSeq = ''
				if ID in protein:
					proteinSeq = protein[ID]
					
					
				### patch the error where cds length is larger than cDNA
				if len(cdsSeq) > len(cdnaSeq):
					cdnaSeq = cdsSeq
				
				blast_hit = ''+'\t'+''
				if ID in homology:
					blast_hit = homology[ID]
					
				
				### add the exon count/length/repeat overlap parameters
				exon = ''+'\t'+''+'\t'+''
				if ID in exons:
					exon = exons[ID]
				
				### add orthoMCL column
				
				oMCL = ' '
				if ID in orthoMCL:
					oMCL = orthoMCL[ID]
				
				### add orthoMCL gene familes
				ofamily = ''+'\t'+''+'\t'+''+'\t'+''+'\t'+''+'\t'+''+'\t'+''
				if oMCL in orthofamily:
					ofamily = orthofamily[oMCL]
				
				### add blast hits
				
				blast_hit_UniProt = ''+'\t'+''+'\t'+''
				if ID in blast:
					blast_hit_UniProt = blast[ID]
				
				### check if protein partial or complete
				if len(proteinSeq) > 2:
					if (proteinSeq[0]=='M') & (proteinSeq[-1]=='*') :
						proteins = 'Complete'
					if (proteinSeq[0]!='M') | (proteinSeq[-1]!='*'):
						proteins = 'Partial'
				if len(proteinSeq) <= 2:
					proteins = 'None'
					
					
				### write the LjRep blast hit
				blast_hit_Ljrep = ''+'\t'+''+'\t'+''
				if ID in ljrblast:
					blast_hit_Ljrep = ljrblast[ID]
				
				print str(count)+'\t'+source+'\t'+ID+'\t'+token[0]+'\t'+token[3]+'\t'+token[4]+'\t'+token[6]+'\t'+cdnaSeq+'\t'+str(len(cdnaSeq))+'\t'+cdsSeq+'\t'+str(len(cdsSeq))+'\t'+proteinSeq+'\t'+str(len(proteinSeq))+'\t'+blast_hit+'\t' + exon + '\t' + str(round(float(len(cdsSeq))/len(cdnaSeq),2)) + '\t' + oMCL + '\t' +ofamily+'\t'+ blast_hit_UniProt + '\t' + proteins +'\t' + gene_name +'\t' + blast_hit_Ljrep


if __name__ == "__main__":
    
    inf, gfile, cfile, pfile, hfile, exon_repeat, oMCLfile, oMCLfamily, blastfile, ljrep_blast = options(sys.argv[1:])
    
    
    ### load cDNA file
    cDNA = hashFasta(gfile)
    
    ### load cds file
    cds = hashFasta(cfile)
    
    ### load protein file
    protein = hashFasta(pfile)
    
    ### hash homology file
    homology = hashHomology(hfile)
    
    ### hash exon_repeat file
    exons = exonRepeat(exon_repeat)
    
    ### hash the orthoMCL
    orthoMCL = hashOrthoMCL(oMCLfile)
    
    ### hash orthoMCL family
    orthofamily = hashOrthofamily(oMCLfamily)
    
    ### hash the blast hits
    blast = hashblast(blastfile)
    
    ### hash the LjRep blast output
    ljrblast = hashljrblast(ljrep_blast)
    
    ### parse file to get the gene names
    format(inf, cDNA, cds, protein, homology, exons, orthoMCL, orthofamily, blast, ljrblast)
