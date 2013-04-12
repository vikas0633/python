### 104_intersect_files_column.py - /Users/vikas0633/Desktop/script/python/ - script to print the desired columns given keys from the files

import os,sys,getopt, re

### main argument to 

def options(argv):
    file1 = ''
    file2 = ''
    key1 = ''
    key2 = ''
    col1 = ''
    col2 = ''
    sep='\t'
    first_line = False
    try:
        opts, args = getopt.getopt(argv,"hfi:j:c:d:s:k:l:",["file1=","file2=","col1=","col2=",'separatedBy=',"key1=","key2="])
    except getopt.GetoptError:
        print '''
            python 104_intersect_files_column.py
                -i <file1>
                -j <file2>
                -c <col1> # multiple columns separated by commas from file 1
                -d <col2> # multiple columns separated by commas from file 2
                -s <separatedBy>
                -k <key1> ### key column from the file 1
                -l <key2> ### key column from the file 2 
                -f <first_line>
                -h <help> 
            
            Example: nice -n 19 python ~/Desktop/script/python/104_intersect_files_column.py \
            --file1 03_gene_association.tair \
            --file2 02_gene_ontology_ext.obo.txt.out \
            --col1  11,5,7,9 \
            --col2  1,2,3,4 \
            --key1 5 \
            --key2 1   
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
            python 104_intersect_files_column.py
                -i <file1>
                -j <file2>
                -c <col1> # multiple columns separated by commas from file 1
                -d <col2> # multiple columns separated by commas from file 2
                -s <separatedBy>
                -k <key1> ### key column from the file 1
                -l <key2> ### key column from the file 2 
                -f <first_line>
                -h <help>
            
        Example: nice -n 19 python ~/Desktop/script/python/104_intersect_files_column.py \
            --file1 03_gene_association.tair \
            --file2 02_gene_ontology_ext.obo.txt.out \
            --col1  11,5,7,9 \
            --col2  1,2,3,4 \
            --key1 5 \
            --key2 1
                '''

            sys.exit(2)
        elif opt in ("-i", "--file1"):
            file1 = arg
        elif opt in ("-j", "--file2"):
            file2 = arg
        elif opt in ("-c", "--col1"):
            col1 = arg
        elif opt in ("-d", "--col2"):
            col2 = arg
        elif opt in ("-k", "--key1"):
            key1 = arg
        elif opt in ("-l", "--key2"):
            key2 = arg
        elif opt in ("-s", "--separatedBy"):
            sep = arg
        elif opt in ("-f","--first_line"):
            first_line = True
    
    return file1, file2, col1, col2, key1, key2, sep, first_line
    

### hash the first file
def HASH(file1,c1,key1,sep,first_line):
    hash = {}
    for line in open(file1,'r'):
        
        if len(line) > 0:
            if (line[0] != '#') & (line[0] != '!'):
                if first_line == False:
                    line = line.strip()
                    token = line.split(sep)
                    lis = list(col1.split(','))
                    value = ''
                    for i in lis:
                        # value += token[int(i)-1]+'\t'
                        value += token[int(i)-1].split('|')[0]+'\t' ### for sake of gene IDs in the gene association file
                    hash[token[int(key1)-1]] = value.strip()
                first_line = False
    return hash


### parse the second file
def PARSE(file2,c2,key2,sep,first_line,hash):
    for line in open(file2,'r'):
        if len(line) > 1:
            if (line[0] != '#') & (line[0] != '!'):
                if first_line == False:
                    line = line.strip()
                    token = line.split(sep)
                    lis = list(col2.split(',')) 
                    value = ''
                    for i in lis:
                        value += token[int(i)-1].split('|')[0]+'\t'
                    if token[int(key2)-1] in hash:
                        print hash[token[int(key2)-1]]+ '\t' + value.strip()
                first_line = False


if __name__ == "__main__":
    
    file1, file2, col1, col2, key1, key2, sep, first_line = options(sys.argv[1:])
    
    ### hash the first file
    hash = HASH(file1,col1,key1,sep,first_line)
    
    ### parse the second file
    PARSE(file2,col2,key2,sep,first_line, hash)
