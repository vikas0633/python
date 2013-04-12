### 12an.py - script to add an identifier based non-redundancy - /Users/vgupta/Desktop/script/python
 
import os,sys,getopt, re

### Usage
'''
python ~/Desktop/script/python/12an.py -f mature_query.txt -i LjMiRCan
'''


### main argument to 

def options(argv):
    file = ''
    identifier = ''
    try:
       opts, args = getopt.getopt(argv,"hf:i:",["file=","identifier="])
    except getopt.GetoptError:
    	print '''
    				python 12an.py 
    					-f <file>
    					-i <identifier> 
    			'''       
    	sys.exit(2)
    
    for opt, arg in opts:
    	if opt == '-h':
    		print '''
    				python 12an.py 
    					-p <file>
    					-s <identifier> 
    			'''
    		sys.exit()
        elif opt in ("-f", "--file"):
        	file = arg
    	elif opt in ("-i", "--identifier"):
        	identifier = arg
        	
    return file, identifier
    
def hash_lines():
	hash = {}
	first_line = True
	for line in open(file,'r'):
		line = line.strip()
		if first_line == False:	
			if line in hash:
				hash[line] += 1
			else:
				hash[line] = 1
		first_line = False
	
	hash_new = {}
	count = 0
	for key in hash:
		count += 1
		iden = identifier +'_'+str(str(count).zfill(3))
		hash_new[key] = iden
	
	return hash,hash_new

def process():
	first_line = True
	for line in open(file,'r'):
		line = line.strip()
		if first_line == False:
			print line + '\t' + str(hash_new[line])+ '\t' + str(hash[line])
		else:
			print 'Sequence'+'\t'+'identifier'+'\t'+'Redundancy'
		first_line = False
		
    
if __name__ == "__main__":
    
    file, identifier = options(sys.argv[1:])
    
    ### add identifier
    hash,hash_new = hash_lines()
    
    ### print identifier
    process()