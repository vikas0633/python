### 21s_summary_eval.py - script for summarizing eval output - /Users/vikas0633/Desktop/script/python
###  python ~/Desktop/script/python/21s_summary_eval.py -i ncbi_321_Lotus.fasta.blastn.gtf_05_transcripts.evm.gtf

import os,sys, getopt
### function to calculate length of each seq in fasta file


def options(argv):
	infile = ''
	try:
		opts, args = getopt.getopt(argv,"hi:",["infile="])
	except getopt.GetoptError:
		print 'python 21s_summary_eval.py -i <infile>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'python 21s_summary_eval.py -i <infile>'
			sys.exit()
		elif opt in ("-i", "--infile"):
			infile = arg
	
	return infile
	
	
def print_summary():
	header = ''
	data = ''
	ln = 0
	for line in open(infile,'r'):
		ln += 1
		line = line.strip()
		token = line.split('\t')
		if ln == 5:
			header += 'Predictions\t'
			data += token[1]+'\t'
		if (7<=ln<=14):
			header += token[0]+'\t'
			data += token[1]+'\t'
			
		if ln == 21:
			header += 'Gene_counts\t'
			data += token[2]+'\t'
		
		if ln == 23:
			header += 'Transcript per Gene\t'
			data += token[2]+'\t'
			
		if ln == 27:
			header += 'Average Length [transript]\t'
			data += token[2]+'\t'
			
		if ln == 29:
			header += 'Total Length [transript]\t'
			data += token[2]+'\t'
		
		if ln == 30:
			header += 'Average coding Length [transript]\t'
			data += token[2]+'\t'
		
		if ln == 32:
			header += 'Total coding Length [transript]\t'
			data += token[2]+'\t'
			
		if ln == 35:
			header += 'Average exon per  [transript]\t'
			data += token[2]+'\t'
		
		if ln == 37:
			header += 'Total exon per  [transript]\t'
			data += token[2]+'\t'
			
		### exon details
		
		if ln == 79:
			header += 'Total count [exon]\t'
			data += token[2]+'\t'
			
		if ln == 80:
			header += 'Average length [exon]\t'
			data += token[2]+'\t'
		
		if ln == 82:
			header += 'Total length [exon]\t'
			data += token[2]+'\t'
		
		### intron details
		if ln == 128:
			header += 'Total Count [intron]\t'
			data += token[2]+'\t'
		
		if ln == 129:
			header += 'Average Length [intron]\t'
			data += token[2]+'\t'
		
		if ln == 131:
			header += 'Total Length [intron]\t'
			data += token[2]+'\t'
		
	
	print header
	print data
		

if __name__ == "__main__":

	infile = options(sys.argv[1:]) 
	
	print_summary()