#27_summary_MS_hit.py - /Users/vikas0633/Desktop/script/python/ - script for process MS hit text file

import sys, re

### open the file
'''
Accession Nr	Sequence	Gene Prediction - Stego_Venom					
		Score	Unique Peptides	Start	End	Duplicates	Peptide Sequences
0_05Tr_v1_fr1	LNVNSHTFIYAIVIFLSVKKSEXNLFLSSLKKNFIRFFLQTKIXQWHKXTCESXHLNKSVFKCDLSXNGFXNQLKLNLMRLALYEFXDDGRSXIMKSFKMLLFPSXSSRGFHPAFLXQLGTRXDPFFSYWKYRILHLHSRCAHHVCLGLGDXILQQCLPKPSCPWRVXEEQCERKLGWWYPSXLGRKSRNPSGHHERMATLIQXYSQLXLDVQRQPSSHRLFAWKXYASDPLHSPTPRRFWHHYGRYHALWANSGRHWRLARTYHLLXRGNDHRXAAAELLCPYQSKDSICLXVHLQMLKIQRKLCLPFLCSXLRSGRGXKHIPSNCDHIAPWXKQQILLWDPRXPKCWYHRNRXEDACQXAXIPHGSLRXWGXKWERRACML	23	1	260	269	1	R.LARTYHLLSR.G
0_05Tr_v1_fr5	EAYTLFAPIFNPIIEDYHEGFKPTDKHPPTDFGDTNTLVNVDPTGEFVVSTRVRCGRSLKGYAFNPCLTEANYKEMEDKVSAVFSAFEGELKGKYYPLTGMDKATQQQLIDDHFLFKEGDRFLQAANACRYWPTGRGIFHNDAKTFLVWVNEEDHLRIISMQKGGDLKAVFERLVKAVNIIESKLPFSRDDRLGFLTFCPTNLGTTIRASVHIALPKLAKDKKVLEDIAAKFNLQVRGTRGEHTESEGGVYDISNKRRMGLTEYQAVKEMQDGILEMIKMEKAASXNFSXFNCAHHPRTHTMQVALNSISIDSKSRFTKGHIXKQIYLNVNSHTFIYAIVIFLSVKKSEXNFFLENLKINFIRFFLQTKIXQWHKXTCESXHL	40	3	122	130	1	R.FLQAANACR.Y
				177	184	1	K.AVNIIESK.L
				259	268	0	R.MGLTEYQAVK.E
'''

proteins = {}
unique_proteins = {}
count_1 = {}
no_x = {}
count = 0
count_2 = {}
count_2_unique = {}
x = {}
withX = {}
for line in open(sys.argv[1],'r'):
	count += 1
	if count > 2:
		token=line.split('\t')
		### load the header name
		if len(token[0]) > 1:
			key = token[0]
			sequence = token[1]
			peptides = int(token[3])
			proteins[token[0]] = ''
			unique_proteins[token[0].split('_')[0]] = ''
			if int(token[3]) > 1: ### check if there are more than 1 peptide hits
				count_1[token[0].split('_')[0]] = ''
				if re.search('X',token[1][int(token[4])-1:int(token[5])+1]):	### make sure x is not present
					x[token[0].split('_')[0]] = ''
				else:
					no_x[token[0].split('_')[0]] = ''
		
		if len(token[4]) > 0:
			if re.search('X',sequence[int(token[4])-1:int(token[5])+1]):	### count if there is a sequence without X in the peptide co-ordinates
				withX[key] = ''
			else:
				count_2_unique[key.split('_')[0]] = ''
				count_2[key] = ''
			
			if peptides > 1: ### unique peptides with at least one sequences containing more than 2 peptides
				if re.search('X',sequence[int(token[4])-1:int(token[5])+1]):	### make sure x is not present
					x[key.split('_')[0]] = ''
				else:
					no_x[key.split('_')[0]] = ''
	
	
print '# proteins hits in the file: '+str(len(proteins))
print '# unique protein hits in the file: '+str(len(unique_proteins))
print '# unique protein hits with more than one peptides: ' + str(len(count_1))
print '# unique protein hits with more than one peptides and with X in at least one: ' + str(len(x))
print '# protein hits with more than one peptides and no X in at least one: ' + str(len(no_x))
print '# proteins hits in the file with at least one peptide containing no X: '+str(len(count_2))
print '# proteins hits in the file with at least one peptide containing X: '+str(len(withX))
print '# unique proteins hits in the file with at least one peptide containing no X: '+str(len(count_2_unique))

