#21v2_format_gff3.py- script to format gff3 file in order to put in MySQL table sequal to 21v - /Users/vikas0633/Desktop/script/python

# Usage: python /Users/vikas0633/Desktop/script/python/21v2_format_gff3.py -i sample.gff3
'''
Usage: python /u/vgupta/script/python/21v2_format_gff3_spider.py \
-i 20120218_Stegodyphous_genemodels.refined.3.mRNA.gff3 \
-g 20130210_Stegodyphous_cdna.refined.fa \
-c 20130210_Stegodyphous_cds.refined.fa \
-p 20130210_Stegodyphous_cds_protein.refined.fa \
-s 20130210_Stegodyphous_cds_protein.refined.fa.blastout.sorted.genewise \
-a 20130210_IPRScan.out.refined \
-e /u/vgupta/02_spider/19_databases4proteomics/Genes_Stego_WB_Augustus_Cufflinks_WBTransDeno.refined.txt \
> 20130218_Stegodyphous_genemodel

																
'''

import os,sys,getopt, re
from A_hash_file import *


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
    gene_densFile = ''
    iprScan_file = ''
    exon_repeat = ''
    oMCLfile = ''
    peptideFile = ''
    try:
       opts, args = getopt.getopt(argv,"hi:g:c:p:b:r:o:f:s:j:d:a:e:", ["ifile=","gfile=","cfile=","pfile=", "bfile=","rfile=","orthoMCL=","oMCLfamily=","blast_search=","ljr_repeat=", "gene_density=", "iprScan=","peptide="])
    except getopt.GetoptError:
       print '''
       			python 21v2_format_gff3.py 
       			-i <inputfile> 
       			-g <gfile> 
       			-c <cfile> 
       			-p <pfile> 
       			-b <homology> 
       			-r <rfile> 
       			-o <orthoMCL> 
       			-f <oMCLfamily> 
       			-s <blast_search> 
       			-j <ljr_repeat> 
       			-d <gene_density> 
       			-a <iprScan> 
       			-e <peptide>
       		''' 
       sys.exit(2)
    for opt, arg in opts:
    	if opt == '-h':
    		print '''
       			python 21v2_format_gff3.py 
       			-i <inputfile> 
       			-g <gfile> 
       			-c <cfile> 
       			-p <pfile> 
       			-b <homology> 
       			-r <rfile> 
       			-o <orthoMCL> 
       			-f <oMCLfamily> 
       			-s <blast_search> 
       			-j <ljr_repeat> 
       			-d <gene_density> 
       			-a <iprScan> 
       			-e <peptide>
       		'''
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
        elif opt in ("-d", "--gene_density"):
        	gene_densFile = arg
        elif opt in ("-a", "--iprScan"):
        	iprScan_file = arg
        elif opt in ("-e", "--peptide"):
        	peptideFile = arg
        
        
    return inputfile, genomefile, cdsfile, proteinfile, homologyfile, exon_repeat, oMCLfile, oMCLfamily, blastfile, ljrep_blast, gene_densFile, iprScan_file, peptideFile

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
	
def hash_denseGene(file):
	''' this function for making a dictionary to store dense genome co-ordinates'''
	dense = {}
	for line in open(file,'r'):
		line = line.strip()
		token = line.split(' ')
		for i in range(int(token[1]), int(token[2])+1,1):
			dense[token[0], i] = token[3]
	return dense


def IPRSCAN(file):
	hash = {}
	for line in open(file,'r'):
		line = line.strip()
		token = line.split('\t')
		hash[token[0]] = token[1]
	return hash

def format(inf, cDNA, cds, protein, homology, exons, orthoMCL, orthofamily, blast, ljrblast, denseGene, iprScan, peptide):
	count = 0
	### print header
	print 'Number'+'\t'+'Source'+'\t'+'ID'+'\t'+'Chromosome'+'\t'+'Start'+'\t'+'End'+'\t'+'Strand'+'\t'+'cDNASeq'+'\t'+'cDNASeq_length'+'\t'+'cdsSeq'+'\t'+'cdsSeq_len'+'\t'+'proteinSeq'+'\t'+'proteinSeq_len'+'\t'+'BlastHit[Medicago]'+'\t'+'BlastHit[Soyabean]'+'\t'+'ExonCount'+'\t'+'ExonicLength'+'\t'+'RepeatFraction' + '\t' + 'CDS/mRNA length ratio' + '\t' + 'OrthoMCL group' + '\t' +'ljr\t'+'mtr\t'+'gmax\t'+'ath\t'+'ptr\t'+'osa\t'+'sob\t'+ 'Annotation [UniProtBlast]' +'\t'+ '%Identity_UniProt' +'\t'+'Evalue_UniProt' + '\t' + 'Protein [Complete/Partial]' + '\tgene_name' + '\t' + 'Annotation [LjRep]' +'\t'+ '%Identity_LjRep' +'\t'+'Evalue_LjRep' + '\t' + 'Gene_density[>5]\t' + 'InterProScan\t' + 'MS-support'
	
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
				elif token[1] == "Stegodyphus_mimosarum.v1.final.fa":
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
					
				### find the gene density
				gene_desnsity = '<5'
				coord = token[0],(int(token[3]) + int(token[4]))/2
				if coord in denseGene:
					gene_desnsity = denseGene[coord]
				
				
				### adding iprScanPro data
				ipr = ' '
				if ID in iprScan:
					ipr = iprScan[ID]
				
				### adding MS hit data
				MS = ''
				tokens = ID.split('.')
				gene_id = '.'.join(tokens[0:len(tokens)-1])
				if gene_id in peptide:
					MS = peptide[gene_id]
				
				
				
				print str(count)+'\t'+source+'\t'+ID+'\t'+token[0]+'\t'+token[3]+'\t'+token[4]+'\t'+token[6]+'\t'+cdnaSeq+'\t'+str(len(cdnaSeq))+'\t'+cdsSeq+'\t'+str(len(cdsSeq))+'\t'+proteinSeq+'\t'+str(len(proteinSeq))+'\t'+blast_hit+'\t' + exon + '\t' + str(round(float(len(cdsSeq))/len(cdnaSeq),2)) + '\t' + oMCL + '\t' +ofamily+'\t'+ blast_hit_UniProt + '\t' + proteins +'\t' + gene_name +'\t' + blast_hit_Ljrep +'\t'+ gene_desnsity + '\t' + ipr + '\t' + MS


if __name__ == "__main__":
    
    inf, gfile, cfile, pfile, hfile, exon_repeat, oMCLfile, oMCLfamily, blastfile, ljrep_blast, gene_densFile, iprScan_file, peptideFile = options(sys.argv[1:])
    
    
    ### load cDNA file
    if gfile != '':
    	cDNA = hashFasta(gfile)
    else:
    	cDNA = {}
    
    ### load cds file
    if cfile != '':
    	cds = hashFasta(cfile)
    else:
    	cds = {}
    
    ### load protein file
    if pfile != '':
    	protein = hashFasta(pfile)
    else:
    	protein = {}
    
    ### hash homology file
    if hfile != '':
    	homology = hashHomology(hfile)
    else:
    	homology = {}
    
    ### hash exon_repeat file
    if exon_repeat != '':
    	exons = exonRepeat(exon_repeat)
    else:
    	exons = {}
    
    ### hash the orthoMCL
    if oMCLfile != '':
    	orthoMCL = hashOrthoMCL(oMCLfile)
    else:
    	orthoMCL = {}
    
    ### hash orthoMCL family
    if oMCLfamily != '':
    	orthofamily = hashOrthofamily(oMCLfamily)
    else:
    	orthofamily = {}
    
    ### hash the blast hits
    if blastfile != '':
    	blast = hashblast(blastfile)
    else:
    	blast = {}
    
    ### hash the LjRep blast output
    if ljrep_blast != '':
    	ljrblast = hashljrblast(ljrep_blast)
    else:
    	ljrblast = {}
    
    
    ### hash the dense gene co-ordinates
    if gene_densFile != '':
    	denseGene = hash_denseGene(gene_densFile)
    else:
    	denseGene = {}
    
    ### hash the IPRscan evidence
    if iprScan_file != '':
    	iprScan = IPRSCAN(iprScan_file)
    else:
    	iprScan = {}
    
    ### hash MS protein support
    if peptideFile != '':
    	peptide = hash_file(peptideFile)
    else:
    	peptide = {}
    	
    ### parse file to get the gene names
    format(inf, cDNA, cds, protein, homology, exons, orthoMCL, orthofamily, blast, ljrblast, denseGene, iprScan, peptide)
