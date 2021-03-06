'''
Created on March 14, 2013

@author: vgupta
'''

#21ai_modify_gene_names.py - script to modify gene names based on N counts - /Users/vikas0633/Desktop/script/python

### check output with 
# python ~/script/python/21ai_modify_gene_names.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame

import os,sys,getopt, re
import A_hash_file
import E_get_chr_size_gff3

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(gff3):
    o.write("Program used: \t\t%s" % "21ai_modify_gene_names.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("GFF3 used: \t\t%s" % gff3+'\n')
            
    

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hg:n:s:",["gff3=","N_counts=","specie="])
    except getopt.GetoptError:
        print '''
            python 21ai_modify_gene_names.py 
                -g <gff3>        [co-ordinate based sorted GFF3]
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
            python 21ai_modify_gene_names.py 
                -g <gff3>        [co-ordinate based sorted GFF3]
                -n <N_counts>    [ File with gap count between genes]
            '''
            sys.exit()
        elif opt in ("-g", "--gff3"):
            gff3 = arg
        elif opt in ("-n", "--N_counts"):
            N_counts = arg
        elif opt in ("-s", "--specie"):
            specie = arg
            
    logfile(gff3)
            
    return gff3, N_counts,specie
    

### hash Ns in from the gap count script
def hash_Ns(ifile):
    hash_count = {}
    for line in open(ifile,'r'):
        line = line.strip()
        token = line.split('\t')
        
        g_id = '.'.join(token[0].split('.')[:-1])
        
        if g_id in hash_count:
            if int(token[1]) > hash_count[g_id] :
                hash_count[g_id] = int(token[1])
        else:
            hash_count[g_id] = int(token[1])
        
    return hash_count
        
### modify gene names
def modifyGeneNames(gff3,chr, hash, specie):
    gene_ID_hash = {}
    count = 0
    for line in open(gff3,'r'):
        line = line.strip()
        if len(line) > 1:
            if line[0] != '#':
                token = line.split('\t')
                if token[0] == chr:         ### make sure that we process line with correct chromosome
                    if token[2] == "gene":
                        gene = line
                        gene_ID = re.search(r'ID=.+',line)
                        gene_ID = gene_ID.group().split(';')[0].replace('ID=','') 
                        #gene_type = re.search(r'Type=.+',line)
                        #gene_type = gene_type.group().split(';')[0].replace('Type=','')
                        if gene_ID not in gene_ID_hash:
                            count += 10
                            if int(hash[gene_ID]) <= 2:
                                count += 0
                            elif 2 < int(hash[gene_ID]) <= 10:
                                count += 10
                            elif 10 < int(hash[gene_ID]) <= 100:
                                count += 100
                            elif 100 < int(hash[gene_ID]) <= 1000:
                                count += 1000
                            elif 1000 < int(hash[gene_ID]) <= 10000:
                                count += 10000
                            elif 10000 < int(hash[gene_ID]) <= 100000:
                                count += 10000
                            elif 100000 < int(hash[gene_ID]):
                                count += 10000
                            gene_ID_hash[gene_ID] = count
                            ### new variable for the keeping the AS cufflinks geneID
                            Count = count
                            flag = True
                        else:
                            flag = False                                            ### flag for new gene
                            Count = gene_ID_hash[gene_ID]                 
                    elif token[2] == "mRNA":
                        match = re.search(r'ID=.+',line)                            ### Get the mRNA ID
                        match = match.group().split(';')[0].replace('ID=','')       
                        variant_no_mRNA = match.split('.')[-1]                           ### Get the mRNA variant number
                        
                            
                        if flag == True:
                            ### print gene 
                            lin = '\t'.join(gene.split('\t')[0:len(gene.split('\t'))-1]) ### all information except ID
                            lin += '\t'+"ID="+token[0].replace('scaffold',specie)+'g'+ str(Count).zfill(5)+';Name='+token[0].replace('scaffold',specie)+'g'+ str(Count).zfill(5)+';'
                            #lin += '\t'+"ID="+chr_name[token[0]]+ str(count).zfill(7)+';Type='+gene_type+';'
                            print lin
                            
                        ### print mRNA
                        mRNA_ID = token[0].replace('scaffold',specie)+'g'+ str(Count).zfill(5)+'.'+variant_no_mRNA
                        lin = '\t'.join(line.split('\t')[0:len(line.split('\t'))-1]) ### all information except ID
                        lin += '\t'+"ID="+mRNA_ID+';Parent='+token[0].replace('scaffold',specie)+'g'+str(Count).zfill(5)+";Name="+mRNA_ID+';'
                        print lin
                    
                    else:
                        ### print rest
                        match = re.search(r'ID=.+',line)
                        match = match.group().split(';')[0].replace('ID=','')       
                        variant_no = match.split('.')[-1].replace('exon','')                           ### Get the mRNA variant number
                        lin = '\t'.join(line.split('\t')[0:len(line.split('\t'))-1]) ### all information except ID
                        lin += '\t'+"ID="+mRNA_ID+'.'+token[2]+'.'+variant_no+';Parent='+mRNA_ID+';'
                        print lin
                    

                            
        

if __name__ == "__main__":
    
    gff3, N_counts, specie = options(sys.argv[1:])
    

    ### hash the N counts from the file
    hash = hash_Ns(N_counts)
    
    ### run it by chromosome
    size = E_get_chr_size_gff3.get_size(gff3)
    for chr in sorted(size):
        ### modify gene names
        modifyGeneNames(gff3,chr, hash, specie)
    
    ### close the logfile
    o.close()