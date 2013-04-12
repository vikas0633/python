### 21t_tau.py - script to add ORF to the gff file - /Users/vikas0633/Desktop/script/python


'''
USAGES:
nice -n 19 python /Users/vgupta/Desktop/script/python/21t_tau.py -f /Users/vgupta/Desktop/temp/Ljchr0-6_pseudomol_20120830.chlo.mito.fa -g sample.20121108_merged_gene_models.gff3

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

def file_empty(file):
    count = sum([1 for line in open(file)])
    if count == 0:
        sys.exit(file+' is empty')	



def options(argv):
	infile = ''
	try:
		opts, args = getopt.getopt(argv,"hf:g:",["genome=","gff="])
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
	
	return infile, gff
	
def extract_seq(infile, header, mRNA_st, mRNA_en):
	os.system('echo ">"'+str(header)+' >temp.fa')
	os.system('nice -n 19 fastacmd -d '+infile+' -p F -s '+header+' -L '+str(mRNA_st)+','+str(mRNA_en) +'| tail +2 >>temp.fa')

def exon(count,mRNA_id,mRNA_st,strand,model_no,mRNA_en):
	exons = {}
	exon_no = 0
	file = './temp'+str(count)+'/temp'+str(count)+'.gff'
	for line in open(file,'r'):
		token = line.strip().split('\t')
		if token[2].startswith('CDS'):
			exon_no += 1
			#Ljchr1_pseudomol_20120830	TAU	CDSf	0	861	.	+	.	temp11 Ljchr1_pseudomol_20120830.1.1;evidence ORF;confidence 3.11627906976744;average_depth 1;
			token=line.split('\t')
			start = int(token[3])
			end = int(token[4])
			ID = 'ID='+mRNA_id+'.'+str(exon_no)+';'+'Parent='+mRNA_id+'.exon;' 
			lin = token[0]+'\t'+token[1]+'\t'+'exon'+'\t'+str(start+mRNA_st)+'\t'+str(end+mRNA_st)+'\t'+token[5]+'\t'+'.'+'\t'+token[7]+'\t'+ID
			o.write(lin+'\n')
			exons[exon_no] = start, end

	return exons

def call_CDS(count,mRNA_id,gene_model,mRNA_st,strand,model_no, exons, mRNA_en):
	file = './temp'+str(count)+'/temp'+str(count)+'.cds.fa'
	for line in open(file,'r'):
		if line.startswith('>'):
			# Ljchr1_pseudomol_20120830	Augustus	mRNA	21237	22777	.	+	.	ID=model.g28684.t1;Parent=gene.g28684.t1
			# >Ljchr1_pseudomol_20120830.1.1:cds:+:2366:3503
			token = line.split('.')
			tokens = line.split(':')
			start = int(tokens[3])
			end = int(tokens[4])
			strand = tokens[2]
			ID = 'ID='+mRNA_id+'.'+str(model_no)+';'+'Parent='+mRNA_id+';'
			frame = 0
			if strand == '+':
				if len(exons)==1:
					five_prime_UTR_st =  exons[1][0] + mRNA_st -1
					five_prime_UTR_en = start - exons[1][0] + mRNA_st -1
					ID = 'ID='+mRNA_id+'.'+str('5_prime_UTR')+';'+'Parent='+mRNA_id+';'
					lin = token[0][1:]+'\t'+'TAU'+'\t'+'five_prime_UTR'+'\t'+str(five_prime_UTR_st)+'\t'+str(five_prime_UTR_en)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
					o.write(lin+'\n')
					CDS_st = start
					CDS_en = end
					ID = 'ID='+mRNA_id+'.'+str('CDS.1')+';'+'Parent='+mRNA_id+';'
					lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(CDS_st + mRNA_st)+'\t'+str(CDS_en + 2+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
					o.write(lin+'\n')
					cds_len = CDS_en - CDS_st + 1
					exon_st = exons[1][0]
					exon_en = exons[1][1]
					last_cds_len = 0
					three_prime_UTR_st = CDS_en + 3+ mRNA_st
				else:
				
					exon_length = 0
					introns = 0
					for i in range(1,len(exons)+1):
						exon_length += exons[i][1]-exons[i][0]
						if i > 1:
							introns += last_exon - exons[i][1] -1
						if exon_length >= start:
							j = i
							break
						last_exon =  exons[i][0]
					five_prime_UTR_st =  exons[1][0] + mRNA_st -1
					five_prime_UTR_en = start - exons[1][0] + mRNA_st -1
					ID = 'ID='+mRNA_id+'.'+str('5_prime_UTR')+';'+'Parent='+mRNA_id+';'
					lin = token[0][1:]+'\t'+'TAU'+'\t'+'five_prime_UTR'+'\t'+str(five_prime_UTR_st)+'\t'+str(five_prime_UTR_en)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
					o.write(lin+'\n')
					CDS_st = start 
					CDS_en = exons[1][1]
					ID = 'ID='+mRNA_id+'.'+str('CDS.1')+';'+'Parent='+mRNA_id+';'
					lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(CDS_st+ mRNA_st)+'\t'+str(CDS_en+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
					o.write(lin+'\n')
					cds_len = CDS_en - CDS_st + 1
					frame = 3- (cds_len) % 3
					for i in range(j+1,len(exons)+1):
						exon_st = exons[i][0]
						exon_en = exons[i][1]
						last_cds_len = cds_len
						cds_len += (exon_en - exon_st) + 1
						ID = 'ID='+mRNA_id+'.'+str('CDS.')+str(i)+';'+'Parent='+mRNA_id+';'
				
						if (start <= exon_st) & (end-start > cds_len):
							lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(exon_st+ mRNA_st)+'\t'+str(exon_en+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
							o.write(lin+'\n')
							frame = 3 - (cds_len)%3
						elif end-start <= cds_len:
							lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(exon_st+ mRNA_st)+'\t'+str(exon_st+(end - start - last_cds_len)+2+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
							o.write(lin+'\n')
							three_prime_UTR_st = exon_st+(end - start - last_cds_len) + 3 + mRNA_st
							break

				three_prime_UTR_en = exon_en + mRNA_st + 1
				ID = 'ID='+mRNA_id+'.'+str('3_prime_UTR')+';'+'Parent='+mRNA_id+';'
				lin = token[0][1:]+'\t'+'TAU'+'\t'+'three_prime_UTR'+'\t'+str(three_prime_UTR_st)+'\t'+str(three_prime_UTR_en)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
				o.write(lin+'\n')
			
			if strand == '-':
				if len(exons)==1:
					five_prime_UTR_st =  exons[len(exons)][1] + mRNA_st
					five_prime_UTR_en =  exons[len(exons)][1] - start  + mRNA_st
					ID = 'ID='+mRNA_id+'.'+str('5_prime_UTR')+';'+'Parent='+mRNA_id+';'
					lin = token[0][1:]+'\t'+'TAU'+'\t'+'five_prime_UTR'+'\t'+str(five_prime_UTR_en)+'\t'+str(five_prime_UTR_st)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
					o.write(lin+'\n')
					ID = 'ID='+mRNA_id+'.'+str('CDS.1')+';'+'Parent='+mRNA_id+';'
	
					CDS_st = exons[len(exons)][1] - start -1
					CDS_en = exons[len(exons)][1] - end
					lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(CDS_en+ mRNA_st)+'\t'+str(CDS_st+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
					o.write(lin+'\n')
					cds_len = CDS_st - CDS_en + 1
					frame = 3- (cds_len) % 3
					last_cds_len = 0
					three_prime_UTR_st = exons[1][1]-(end -start -last_cds_len) - 1 + mRNA_st
				else:
					exon_length = 0
					introns = 0
					for i in range(len(exons),0,-1):
						exon_length += exons[i][1]-exons[i][0]
						if i < len(exons):
							introns += last_exon - exons[i][1] -1
						if exon_length >= start:
							j = i
							break
						last_exon =  exons[i][0]
					
					five_prime_UTR_st =  exons[len(exons)][1] + mRNA_st
					five_prime_UTR_en =  exons[len(exons)][1] - start - introns + mRNA_st
					ID = 'ID='+mRNA_id+'.'+str('5_prime_UTR')+';'+'Parent='+mRNA_id+';'
					lin = token[0][1:]+'\t'+'TAU'+'\t'+'five_prime_UTR'+'\t'+str(five_prime_UTR_en)+'\t'+str(five_prime_UTR_st)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
					o.write(lin+'\n')
					ID = 'ID='+mRNA_id+'.'+str('CDS.1')+';'+'Parent='+mRNA_id+';'

							
					CDS_st = exons[len(exons)][1] - start - introns -1
					CDS_en = exons[j][0] 
					lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(CDS_en+ mRNA_st)+'\t'+str(CDS_st+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
					o.write(lin+'\n')
					cds_len = CDS_st - CDS_en + 1
					frame = 3- (cds_len) % 3
					for i in range(j-1,0,-1):
						exon_st = exons[i][0]
						exon_en = exons[i][1]
						last_cds_len = cds_len
						cds_len += (exon_en - exon_st) + 1
						ID = 'ID='+mRNA_id+'.'+str('CDS.')+str(i)+';'+'Parent='+mRNA_id+';'
				
						if (start <= exon_st) & (end-start > cds_len):
							lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(exon_st+ mRNA_st)+'\t'+str(exon_en+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
							o.write(lin+'\n')
							frame = 3- (cds_len) % 3
						elif end-start <= cds_len:
							lin = token[0][1:]+'\t'+'TAU'+'\t'+'CDS'+'\t'+str(exons[i][1]-(end -start -last_cds_len)+mRNA_st)+'\t'+str(exons[i][1]+ mRNA_st)+'\t'+'.'+'\t'+strand+'\t'+str(frame)+'\t'+ID
							o.write(lin+'\n')
							three_prime_UTR_st = exons[i][1]-(end -start -last_cds_len) - 1 + mRNA_st
							break
				three_prime_UTR_en = exons[1][0] + mRNA_st - 1
				ID = 'ID='+mRNA_id+'.'+str('3_prime_UTR')+';'+'Parent='+mRNA_id+';'
				lin = token[0][1:]+'\t'+'TAU'+'\t'+'three_prime_UTR'+'\t'+str(three_prime_UTR_en)+'\t'+str(three_prime_UTR_st)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID
				o.write(lin+'\n')
			
			
			os.system('echo ">" '+ID+' >> Ljr_cds.fa')
			os.system('cat '+'./temp'+str(count)+'/temp'+str(count)+'.cds.fa | tail +2 >> Ljr_cds.fa')
			os.system('echo ">" '+ID+' >> Ljr_cds_protein.fa')
			os.system('cat '+'./temp'+str(count)+'/temp'+str(count)+'.protein.fa| tail +2 >> Ljr_cds_protein.fa')

			
					
			'''	
			#lin = token[0][1:]+'\t'+gene_model+'\t'+'mRNA'+'\t'+str(start)+'\t'+str(end)+'\t'+'.'+'\t'+strand+'\t'+'.'+'\t'+ID		
			#o.write( lin
			'''

def mRNA(count,gene_id,gene_model,mRNA_st,strand,model_no):
	file = './temp'+str(count)+'/temp'+str(count)+'.cdna.fa'
	for line in open(file,'r'):
		if line.startswith('>'):
			# Ljchr1_pseudomol_20120830	Augustus	mRNA	21237	22777	.	+	.	ID=model.g28684.t1;Parent=gene.g28684.t1
			# Ljchr1_pseudomol_20120830.1.1:cdna:1:1541:+
			token = line.split('.')
			tokens = line.split(':')
			start = int(tokens[2]) + mRNA_st
			end = int(tokens[3]) + mRNA_st
			ID = 'ID='+gene_id+'.'+str(model_no)+';'+'Parent='+gene_id+';'
			mRNA_id = gene_id+'.'+str(model_no)
			lin = token[0][1:]+'\t'+gene_model+'\t'+'mRNA'+'\t'+str(start)+'\t'+str(end)+'\t'+'.'+'\t'+tokens[4].strip()+'\t'+'.'+'\t'+ID		
			o.write(lin+'\n')
			print file, model_no,"flag", lin
			### call exons from TAU gff file
			exons = exon(count,mRNA_id,mRNA_st,strand,model_no,end)
			### call CDS/UTRs
			call_CDS(count,mRNA_id,gene_model,mRNA_st,strand,model_no,exons,end)
	

def process_transcript(count,gene_id,gene_model,mRNA_st,strand,model_no):
	### open cdna file and o.write( as mRNA gene-model
	mRNA(count,gene_id,gene_model,mRNA_st,strand,model_no)
	

def process_gff(gff,infile):
	id = ''
	count = -1
	gene_flag = ''
	new_transcript = False
	for line in open(gff,'r'):
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
			# Ljchr1_pseudomol_20120830	Augustus	mRNA	21237	22777	.	+	.	ID=model.g28684.t1;Parent=gene.g28684.t1
			if new_transcript == True:
				tau_in.close()
				### make genome fasta file
				extract_seq(infile, header, mRNA_st, mRNA_en)
				### change the format of the gff3 file
				os.system('nice -n 19 python /Users/vikas0633/Desktop/script/python/21u_make_gff2.py -i tau_in -o temp.gff')
			
				### run TAU
				os.system('nice -n 19 perl /Users/vikas0633/Desktop/tools/TAU/TAU.pl -A temp.fa -G temp.gff -O temp'+str(count) +' ')	
				
				file = './temp'+str(count)+'/temp'+str(count)+'.cdna.fa'

				if sum([1 for line1 in open(file)]) != 0:
					### process transcript
					process_transcript(count,last_gene_id,last_gene_model,mRNA_st,strand,last_model_no)				
				
				#sys.exit(0)
			new_transcript = True
			header = token[0]
			mRNA_st = int(token[3]) 
			mRNA_en = int(token[4]) 
			strand = token[6]
			start = int(token[3]) - mRNA_st
			end = int(token[4]) - mRNA_st
			lin = token[0]+'\t'+token[1]+'\t'+token[2]+'\t'+str(start)+'\t'+str(end)+'\t'+token[5]+'\t'+token[6]+'\t'+token[7]+'\t'+token[8]
			tau_in = open('tau_in','w')
			id = line.split('=')[1].split(';')[0]
			tau_in.write(lin+'\n')
			if model_no == 1:
				last_gene_id = gene_id
				last_gene_model = gene_model
			last_model_no = model_no
			model_no += 1
			
			
		elif line.split('=')[2].split(';')[0] == id:
			start = int(token[3]) - mRNA_st
			end = int(token[4]) - mRNA_st
			lin = token[0]+'\t'+token[1]+'\t'+token[2]+'\t'+str(start)+'\t'+str(end)+'\t'+token[5]+'\t'+token[6]+'\t'+token[7]+'\t'+token[8]
			tau_in.write(lin+'\n')
		
		
	tau_in.close()
	count += 1
	model_no += 1
	if model_no == 1:
		gene_name = gene_name.replace('+','.').replace('-','.')
		o.write(gene_name+'\n')
	### make genome fasta file
	extract_seq(infile, header, mRNA_st, mRNA_en)
	### change the format of the gff3 file
	os.system('nice -n 19 python /Users/vikas0633/Desktop/script/python/21u_make_gff2.py -i tau_in -o temp.gff')

	### run TAU
	os.system('nice -n 19 perl /Users/vikas0633/Desktop/tools/TAU/TAU.pl -A temp.fa -G temp.gff -O temp'+str(count) +' ')	
	
	file = './temp'+str(count)+'/temp'+str(count)+'.cdna.fa'
	if sum([1 for line in open(file)]) != 0:
		### process transcript
		process_transcript(count,last_gene_id,last_gene_model,mRNA_st,strand,model_no)
		
				
if __name__ == "__main__":

	infile,gff = options(sys.argv[1:]) 
	
	
	### sort the gff file
	os.system('sort -s -nk 4,4 '+gff+" |sort -s -k 1,1|perl -pe '$_=~s/\+/./'|perl -pe '$_=~s/\-/./'| uniq > sorted.gff3")
	
	#file_empty(infile)
	#file_empty(gff)
	
	gff = 'sorted.gff3'
	
	o = open('20121112_TAU_genemodel.gff3','w')
	### process gene model file
	process_gff(gff,infile)
	
	#os.system('rm -r temp*')