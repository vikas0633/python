import sys

'''
Usages: python ~/script/python/26_summary_mirDeepP.py miRNA_predictions miRNA_filter_P_prediction
'''

### hash the miRDeep output file
	
'''
score_star	-1.3
score_randfold	1.6
score_mfe	2.5
score_freq	-1.8
score	1.2
flank_first_end  	194
flank_first_seq  	CAACGCTCAACCCTAATCCTTGGAGTATATAAAGAGGCTGAAGACTGAAAGAAACGGCTAGAAGCATTCACATACGCGCAAGACAAACTTAAATTCTTCTAAGCTTTCTTTCATCTGAAATTCATTGAGTTTACTATTAGTCTTTTAGAAGCAAATCTCTTGTAAACAATTCTTTGATAAACAGTTTGTTTAGT
flank_first_struct  	...(((((.......(((...(((((.....((((((((((.((.((((((((..(((((((((.(((........(....).........))).))))).))))))))))))))..((((((...))))))....))))))))))....((((.....)))).....)))))..)))..........((((((
flank_second_beg  	234
flank_second_seq  	CTAAGAGAGTGAATCTT
flank_second_struct  	))))).)))))......
freq  	8
loop_beg  	211
loop_end  	214
loop_seq  	GGCT
loop_struct  	...)
mature_arm  	second
mature_beg  	215
mature_end  	233
mature_query  	GATCGGATCCTAGAGAAGA_x5
mature_seq  	GATCGGATCCTAGAGAAGA
mature_strand  	+
mature_struct  	))))...))))))))..))
pre_seq  	TCCTTTAGGAGATCAAGGCTGATCGGATCCTAGAGAAGA
pre_struct  	(.(((((((((((((....)))))...))))))))..))
pri_beg  	1
pri_end  	250
pri_id  	chr0_10484
pri_mfe  	-59.83
pri_seq  	CAACGCTCAACCCTAATCCTTGGAGTATATAAAGAGGCTGAAGACTGAAAGAAACGGCTAGAAGCATTCACATACGCGCAAGACAAACTTAAATTCTTCTAAGCTTTCTTTCATCTGAAATTCATTGAGTTTACTATTAGTCTTTTAGAAGCAAATCTCTTGTAAACAATTCTTTGATAAACAGTTTGTTTAGTTCCTTTAGGAGATCAAGGCTGATCGGATCCTAGAGAAGACTAAGAGAGTGAATCTT
pri_struct  	...(((((.......(((...(((((.....((((((((((.((.((((((((..(((((((((.(((........(....).........))).))))).))))))))))))))..((((((...))))))....))))))))))....((((.....)))).....)))))..)))..........(((((((.(((((((((((((....)))))...))))))))..))))))).)))))......
star_arm  	first
star_beg  	195
star_end  	210
star_seq  	TCCTTTAGGAGATCAA
star_struct  	(.(((((((((((((.
TGATCGGATCCTAGAGAAG_x3	19	1..19	chr0_10484	250	214..232	1e-04	1.00	42.1	Plus / Plus
GATCGGATCCTAGAGAAGA_x5	19	1..19	chr0_10484	250	215..233	1e-04	1.00	42.1	Plus / Plus
'''

def hash(hash_miRDeep):
	hash_miRDeep[pri_id,mature_query] = pri_id+'\t'+mature_query+'\t'+star_seq+'\t'+pri_seq+'\t'+pri_struct

	return hash_miRDeep
	

hash_miRDeep = {}
first_line = True
for line in open(sys.argv[1],'r'):
	line = line.strip()
	token = line.split()
	if line.startswith('mature_query'):
		mature_query = token[1]
	if line.startswith('pri_id'):
		pri_id = token[1]
	if line.startswith('pri_seq'):
		pri_seq = token[1]
	if line.startswith('pri_struct'):
		pri_struct = token[1]
	if line.startswith('star_seq'):
		star_seq = token[1]
	if first_line == False:
		if line.startswith('score_star'):
			hash_miRDeep = hash(hash_miRDeep)
	first_line = False

hash_miRDeep = hash(hash_miRDeep)	

### Process the final miRNA prediction output
'''
chr1	+	GCTTATAAATAGGACCGGA_x38	chr1_10389	19643276..19643294	19643276..19643401	GCTTATAAATAGGACCGGA	GCTTATAAATAGGACCGGAGGGAGTACCATATATGTATGTCTCGGTGAGCGGTGTAAAGAAATAATTGTTCTCTTTATATTCATGTGTGTTTCTTAAATACTCCATCCGTTCCTATTTATAAGAAC
'''

print 'Pri_id'+'\t'+'mature_query'+'\t'+'star_seq'+'\t'+'pri_seq'+'\t'+'pri_struct'+'\t'+'Chromosome'+'\t'+'Start..End'

for line in open(sys.argv[2],'r'):
	line = line.strip()
	token = line.split()
	
	print hash_miRDeep[token[3],token[2]]+'\t'+token[0]+'\t'+token[4]