### E_get_chr_size_gff3.py - script takes a gff3 file and returns max position for each chromosome - /Users/vikas0633/Desktop/script/python/

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
    return size