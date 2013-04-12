### 20d_count_mapped_fastq_inSam.py - script for counting common reads in two fastq files - script for counting the reads mapped



import os,sys, getopt

def file_empty(file):
    count = sum([1 for line in open(file)])
    if count == 0:
        sys.exit(file+' is empty')	



def options(argv):
	infile = ''
	try:
		opts, args = getopt.getopt(argv,"hf:s:",["fastq=","sam="])
	except getopt.GetoptError:
		print 'python  20d_count_mapped_fastq_inSam.py -f <fastq> -s <sam>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'python  20d_count_mapped_fastq_inSam.py -f <fastq> -s <sam>'
			sys.exit()
		elif opt in ("-f", "--fastq"):
			fastq = arg
		elif opt in ("-s", "--sam"):
			sam = arg
	
	return fastq, sam

### hash the sam file
def hash_sam(sam):
	ln = 0
	hash = {}
	for line in open(sam,'r'):
		ln += 1
		token = line.strip().split('\t')
		hash[token[0]] = ''
	
		if ln%10000 == 0:
			print ln,"Sam lines processed"
	
	### sam file hashed
	return hash

### hash the fastq files
def process_fastq(fastq,hash):
	count = 0
	ln = 0
	for line in open(fastq,'r'):
		ln += 1
		if ln%4 == 2:
			line = line.strip()
			key1 = line.split(' ')[0]
			key2 = line.split('/')[0]
			if key1 in hash:
				count += 1
			elif key2 in hash:
			 	count += 1
		if ln%100000 == 0:
			print ln,"fastq lines processed"		 	
	print "fastq read count: ", ln
	print "Mapped read count: ", count
	print "Fraction of reads mapped", float(count)*100/ln
		


if __name__ == "__main__":
	
	fastq, sam = options(sys.argv[1:]) 
	
	### check if files are empty
	'''
	file_empty(fastq)
	file_empty(sam)
	'''
	### hash the sam file
	hash = hash_sam(sam)
	
	### process the fastq file
	
	
	fastq="/home/vgupta/01_genome_annotation/02_transcriptomics_data/2010_02_17_Fasteris_MG20_Gifu_transcripts/100128_s_1_1_seq_GHD-1.txt"
	print fastq
	process_fastq(fastq,hash)
	
	fastq="/home/vgupta/01_genome_annotation/02_transcriptomics_data/2010_02_17_Fasteris_MG20_Gifu_transcripts/100128_s_2_1_seq_GHD-2.txt"
	print fastq
	process_fastq(fastq,hash)
	
	fastq="/home/vgupta/01_genome_annotation/02_transcriptomics_data/2010_03_22_Fasteris_MG20_Gifu_transcripts/100226_s_7_1_seq_GHD-1.txt"
	print fastq
	process_fastq(fastq,hash)
	
	fastq="/home/vgupta/01_genome_annotation/02_transcriptomics_data/2010_03_22_Fasteris_MG20_Gifu_transcripts/100226_s_8_1_seq_GHD-2.txt"
	print fastq
	process_fastq(fastq,hash)
	
	fastq="/home/vgupta/01_genome_annotation/02_transcriptomics_data/2012_05_16_kama_mRNAseq/Sample_P0036_N034-01/P0036_N034-01_CAGATC_L008_R1_001.fastq"
	print fastq
	process_fastq(fastq,hash)
		
	fastq="/home/vgupta/01_genome_annotation/02_transcriptomics_data/2012_05_16_kama_mRNAseq/Sample_P0036_N034-02/P0036_N034-02_ACTTGA_L008_R1_001.fastq"
	print fastq
	process_fastq(fastq,hash)
		
	fastq="/home/vgupta/01_genome_annotation/02_transcriptomics_data/2012_05_16_kama_mRNAseq/Sample_P0036_N034-03/P0036_N034-03_GATCAG_L008_R1_001.fastq"
	print fastq
	process_fastq(fastq,hash)
	
	fastq="/home/vgupta/01_genome_annotation/02_transcriptomics_data/2012_05_16_kama_mRNAseq/Sample_P0036_N034-04/P0036_N034-04_TAGCTT_L008_R1_001.fastq"
	print fastq
	process_fastq(fastq,hash)	
	
	fastq="/home/vgupta/01_genome_annotation/02_transcriptomics_data/2012_05_16_kama_mRNAseq/Sample_P0036_N034-05/P0036_N034-05_GGCTAC_L008_R1_001.fastq"
	print fastq
	process_fastq(fastq,hash)
		
	fastq="/home/vgupta/01_genome_annotation/02_transcriptomics_data/2012_05_16_kama_mRNAseq/Sample_P0036_N034-06/P0036_N034-06_CTTGTA_L008_R1_001.fastq"
	print fastq
	process_fastq(fastq,hash)	
	
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
	
	
	
	
	
	
	