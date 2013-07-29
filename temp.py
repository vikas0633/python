import sys, re

### get ID
def get_ID(line):
    line = line.strip()
    match = re.search(r'ID=.+;',line)
    if match:
        return match.group().split(';')[0].replace('ID=','')
    
### get ID
def get_PARENT(line):
    line = line.strip()
    match = re.search(r'Parent=.+;',line)
    if match:
        return match.group().split(';')[0].replace('Parent=','')
    
gene_names = {}
for line in open(sys.argv[1],'r'):
    line = line.strip()
    if len(line) > 1 :
        if not line.startswith('#'):
            token = line.split('\t')
            if token[2] == "gene":
                g_id = get_ID(line)
                gene_names[g_id] = ''
            
for line in open(sys.argv[1], 'r'):
    line = line.strip()
    if len(line) > 1 :
        if not line.startswith('#'):
            token = line.split('\t')
            if token[2] == "mRNA":
                g_id = get_PARENT(line)
                if g_id not in gene_names:
                    print 'Parent Missing'
                    print line
                    sys.exit()