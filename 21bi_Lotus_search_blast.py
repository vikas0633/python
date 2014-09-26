#-----------------------------------------------------------+
#                                                           |
# template.py - Python template                             |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                        |
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
import os,sys,getopt, re

import threading

from multiprocessing import Process, Queue, Manager
from threading import Thread
import classGene
### global variables
global infile

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')


def get_size(file):
    size = {}
    for line in open(file,'r'):
        line = line.strip()
        if len(line) > 0:
            if line[0] != '#':
                token = line.split('\t')
                
                if token[0] not in size:
                    size[token[0]] = int(token[4])
                    
                else:
                    if int(token[4]) > size[token[0]]:
                        size[token[0]] = int(token[4])
    size_sorted={}
    for w in sorted(size, key=size.get, reverse=False):
        size_sorted[w]=size[w]
    return size_sorted


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

def temp(file):
    return

def LOADfasta(file):
	first_line = True
	seq = {}
	string = ''
	for line in open(file,'r'):
		line = line.strip()
		if len(line) > 0 :			
			if line[0] == '>':
				if first_line == False:
					if string != '': 
						seq[header] = string
				string = ''
				header = '.'.join(line[1:].strip().split()[0].split('.')[:-1])
			else:
				string += line
		first_line = False			
	if string != '': 
		seq[header] = string
	return seq

def hash_anno():
    AnnoHash = {}
    for line in open(anno, 'r'):
        line = line.strip()
        tokens = line.split('\t')
        g_id = '.'.join(tokens[0].split('.')[:-1])
        AnnoHash[g_id] = tokens[1].replace(' ','_')
    
    return AnnoHash

def create_blast_database(database):
    os.system('nice -n 19 makeblastdb -in '+ database +" -dbtype 'prot'")

### main argument to
def options(argv):
	global infile, threads, anno, DATABASE1, DATABASE2
	infile = ''
	threads = 2
	try:
		opts, args = getopt.getopt(argv,"hi:t:a:d:e:",["infile=","threads=","anno=","DATABASE1=","DATABASE2="])
	except getopt.GetoptError:
		help()
	for opt, arg in opts:
		if opt == '-h':
			help()
		elif opt in ("-i", "--infile"):
			infile = arg
		elif opt in ("-t", "--threads"):
			threads = int(arg)
		elif opt in ("-a", "--anno"):
			anno  = arg
		elif opt in ("-d", "--DATABASE1"):
			DATABASE1  = arg
		elif opt in ("-e", "--DATABASE2"):
			DATABASE2  = arg
	
	logfile(infile)

def blast_seq(name, seq, DATABASE):
    ### create temp file with sequence
    o = open('temp.fa','w')
    o.write('>'+name+'\n')
    o.write(seq+'\n')
    o.close()
    os.system('nice -n 19 blastp -num_threads '+ str(threads) +' -max_target_seqs 500 -db '+ DATABASE+ ' -query temp.fa -outfmt 6 -out temp.out')


def process_seq(header, seq):
	global temp_line
	target_seq = 'No_Hit'
	blast_seq(header, seq, DATABASE2) ## Lj -> At -1 
	### extract the best hit
	for line in open('temp.out', 'r'):
		line = line.strip()
		target_seq = '.'.join(line.split('\t')[1].split()[0].split('.')[:-1]) #At top hit
		break
    
	temp_line += '\t' + target_seq 
	
	if target_seq != 'No_Hit':
		added = False ### flag for adding position
		sec_target_seq = ''
		### fetch target seq from Database 2
		first_line = True
		seq = HASH_DATABASE2[target_seq]
		last_target_seq = ''
		### blast against the first database ## At -> Lj -2
		blast_seq(target_seq, seq, DATABASE1)
		At_header = target_seq
		count = 0
		### find the position of original query
		for line in open('temp.out', 'r'):
			line = line.strip()
			if first_line == True:
				sec_target_seq = '.'.join(line.split('\t')[1].split()[0].split('.')[:-1]) #Lj top hit
				print >> sys.stderr, "Found second target: ", sec_target_seq
			target_seq = line.split('\t')[1].split()[0] ## Lj pos
			if len(line) > 0 and last_target_seq != target_seq:
				count += 1
			if re.search(header, target_seq):
				temp_line += '\t'+ str(count)
				added = True
				break
			last_target_seq = target_seq
			first_line = False
		if added == False:
			temp_line += '\t'+ str('No_Hit')
		
		temp_line += '\t' + sec_target_seq 
		
		if sec_target_seq != 'No_Hit':
			added = False ### flag for adding position
			### fetch target seq from Database 1
			seq = HASH_DATABASE1[sec_target_seq]
			last_target_seq = ''
			### blast against the second database
			blast_seq(sec_target_seq, seq, DATABASE2) ## Lj -> At -3
			count = 0
			### find the position of original query
			for line in open('temp.out', 'r'):
				line = line.strip()
				target_seq = line.split('\t')[1].split()[0] ## At pos
				if len(line) > 0 and last_target_seq != target_seq:
					count += 1
				if re.search(At_header, target_seq):
					temp_line += '\t'+ str(count)
					added = True
					break
				last_target_seq = target_seq
			
			if added == False:
				temp_line += '\t'+ str('No_Hit')
			
	return temp_line
	
def read_query():
    global temp_line
    first_line = True
    ### read in the query file and blast individually
    for line in open(infile,'r'):
        line = line.strip()
        
        if len(line) > 0:
            if first_line == False:
				header = re.findall(r'Lj.g.......',line.split()[0])[0] ###  make sure input file agree Lotus pattern
				print >> sys.stderr, "Running gene: " + str(header)
				seq=HASH_DATABASE1[header]
				temp_line=line.split('\t')[0] + '\t'+ header
				print process_seq(header, seq)
        first_line = False
        

if __name__ == "__main__":
	
	options(sys.argv[1:])
	
	start_time = datetime.datetime.now()
	print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())

	### hash annotations
	global AnnoHash, HASH_DATABASE1, HASH_DATABASE2
	AnnoHash = hash_anno()
	
	### create blast databases for two genome
	#create_blast_database(DATABASE1)
	### hash fasta file
	HASH_DATABASE1 = LOADfasta(DATABASE1)
	print >> sys.stderr, "Number of sequences in first database: " + str(len(HASH_DATABASE1))
	#create_blast_database(DATABASE2)
	HASH_DATABASE2 = LOADfasta(DATABASE2)
	print >> sys.stderr, "Number of sequences in secon database: " + str(len(HASH_DATABASE2))
	### blast query sequences
	read_query()
	
	#print >> sys.stderr, "Output count: " + str(temp(infile))
	print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
	print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)
    