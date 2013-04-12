### for sRNA sequences
# script for calculating number of mistaches in Blast hits against PMRD
## updated on 20110101
## uses: 20110506_count_mismatch.py /Users/vikasgupta/Desktop/tools/ncbi-blast-2.2.24+/db/PMRD_all.fa /Users/vikasgupta/Desktop/plant/2011_week4/20110128_PMRD_mapping/PMRD_species_wise_dataset.txt 2010-03-26_GHD-5-16_19-38_all.profile_ext.cut20.norm.scores.fasta_blastout >2010-03-26_GHD-5-16_19-38_all.profile_ext.cut20.norm.scores.fasta_blastout_mismatch_count 

import sys



PMRD = open(sys.argv[1],'r')

#### just for using name of miRNA, nothing to do with mapping
PMRD_seq = {}
for line in PMRD.readlines():
        line = line.strip()
        token=line.split()
        if(line[0]=='>'):
                key = token[0][1:]
        else:
                PMRD_seq[key] = line
                #print(key,(PMRD_seq[key]))
#print len(PMRD_seq)

### inout for couting mismatches
f = open(sys.argv[2],'r')
last_seq = '1'
count_0 = 0
count_1 = 0
count_2 = 0
count_3 = 0
number = 1
### best hit equals best match
score = 0
best_match = "none"
best_family = "none"
for line in f.readlines():
        line =line.strip()
        token = line.split()
        token_part = token[0].split('_')
        #print(line)
        parts = last_seq.split('_')
        if(number < int(parts[0])):
                for i in range(number,int(parts[0])):
                        print (str(i)+'\t'+'0:0:0:0'+'\t'+"none"+'\t'+"none")
                        m = 3
                number = int(parts[0])
        #print(line)
        min_len = min(int(len(PMRD_seq[token[1]])),int(token_part[2]))
        #print(min_len,int(token[3]))
        
        if( last_seq == token[0]):
                ## Case1: No mis-matches
                if( (int(token[3])==int(min_len))&(int(token[4])==0)):
                        #print("match",min_len,int(token[3]),int(PMRD_seq[token[1]]),int(token_part[2]))
                        count_0 += 1
                        score = 1000
                        best_match = PMRD_seq[token[1]]
                        best_family = token[1]
                 ## Case2: Exactly one mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                if((int(token[3])+1-int(token[4]))==int(min_len)):
                        count_1 +=1
                        if ( score < 100):
                                best_match = PMRD_seq[token[1]]
                                best_family = token[1]
                                score = 100
                                #print(best_match,best_family)
                 ## Case3: Exactly two mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                if((int(token[3])+2-int(token[4]))==int(min_len)):
                        count_2 +=1
                        if ( score < 10):
                                best_match = PMRD_seq[token[1]]
                                best_family = token[1]
                                score = 10
               ## Case4: Exactly three mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                if((int(token[3])+3-int(token[4]))==int(min_len)):
                        count_3 +=1
                        if ( score < 1):
                                best_match = PMRD_seq[token[1]]
                                best_family = token[1]
                                score = 1
                else:
                        continue
        else:
                #print score
                print(last_seq+'\t'+str(count_0)+':'+str(count_1)+':'+str(count_2)+':'+str(count_3)+'\t'+best_match+'\t'+best_family)
                score = 0
                count_0 = 0
                count_1 = 0
                count_2 = 0
                count_3 = 0
                best_match = "none"
                best_family = "none"
                method_used = "none"
                number += 1
        ## Case2: One mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                
        
        last_seq = token[0]


"""

### for mirDeep header
# script for calculating number of mistaches in Blast hits against PMRD
total_match_count = 0
PMRD = open('/Users/vikasgupta/Desktop/tools/ncbi-blast-2.2.24+/db/PMRD_all.fa','r')
PMRD_seq = {}
for line in PMRD.readlines():
        line = line.strip()
        if(line[0]=='>'):
                key = line[1:]
        else:
                PMRD_seq[key] = line
                #print(key,(PMRD_seq[key]))

f = open('/Users/vikasgupta/Desktop/plant/2011_week1/20110103/miRNA_list_predictions_PMRD_v5_blast_output.fasta','r')
last_seq = '1'
count_0 = 0
count_1 = 0
count_2 = 0
count_3 = 0
number = 1
### best hit equals best match
score = 0
best_match = "none"
best_family = "none"
for line in f.readlines():
        line =line.strip()
        token = line.split()
        token_part = token[0].split('size')
        parts = last_seq.split('size')
        parts[0]=parts[0].replace('.','')
        if(number < int(parts[0])):
                for i in range(number,int(parts[0])):
                        print (str(i)+'\t'+'0:0:0:0'+'\t'+"none"+'\t'+"none")
                number = int(parts[0])
        #print(line)
        min_len = min(int(len(PMRD_seq[token[1]])),int(token_part[1][0:2]))
        #print(min_len,int(token[3]))
        
        if( last_seq == token[0]):
                ## Case1: No mis-matches
                if( int(token[3])==int(min_len)&(int(token[4])==0)):
                        #print("match",min_len,int(token[3]),int(PMRD_seq[token[1]]),int(token_part[2]))
                        count_0 += 1
                        score = 1000
                        best_match = PMRD_seq[token[1]]
                        best_family = token[1]
                 ## Case2: Exactly one mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                if((int(token[3])+1-int(token[4]))==int(min_len)):
                        count_1 +=1
                        if ( score < 100):
                                best_match = PMRD_seq[token[1]]
                                best_family = token[1]
                                score = 100
                 ## Case3: Exactly two mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                if((int(token[3])+2-int(token[4]))==int(min_len)):
                        count_2 +=1
                        if ( score < 10):
                                best_match = PMRD_seq[token[1]]
                                best_family = token[1]
                                score = 10
                 ## Case4: Exactly three mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                if((int(token[3])+3-int(token[4]))==int(min_len)):
                        count_3 +=1
                        if ( score < 1):
                                best_match = PMRD_seq[token[1]]
                                best_family = token[1]
                                score = 1
                else:
                        continue
        else:
                print(last_seq+'\t'+str(count_0)+':'+str(count_1)+':'+str(count_2)+':'+str(count_3)+'\t'+best_match+'\t'+best_family)
                if( best_match != "none"):
                        total_match_count += 1
                count_0 = 0
                count_1 = 0
                count_2 = 0
                count_3 = 0
                best_match = "none"
                best_family = "none"
                number += 1
        ## Case2: One mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                
        
        last_seq = token[0]

print ( "total_match_count: ", total_match_count)
"""
"""
## for Browsed PMRD data files 
## Update on 20110104
PMRD = open('/Users/vikasgupta/Desktop/tools/ncbi-blast-2.2.24+/db/20110104_PMRD_browseBased_9229.fa','r')
PMRD_seq = {}
PMRD_family = {}
PMRD_method = {}
PMRD_reference = {}
for line in PMRD.readlines():
        line = line.strip()
        if(line[0]=='>'):
                key = line[1:]
        else:
                PMRD_seq[key] = line
a = 0
PMRD_data = open('/Users/vikasgupta/Desktop/plant/2011_week1/20110104/PMRD_species_wise_dataset.txt','r')
for line in PMRD_data.readlines():
        line = line.strip()
        token = line.split('\t')
        if(token[0]=='1'):
                reference = token[5]
        if token[2] in PMRD_seq:
                a += 1
                PMRD_family[token[2]]= token[1]
                PMRD_method[token[2]]= token[4]
                PMRD_reference[token[2]]= reference
                #print(PMRD_family[token[2]],token[2],PMRD_seq[token[2]])

f = open('/Users/vikasgupta/Desktop/plant/2011_week1/20110104/2010_12_29_einf_var_nonred_profiles_file3_blast_output_20110104_Browsed_PMRD.fasta','r')
last_seq = '1'
count_0 = 0
count_1 = 0
count_2 = 0
count_3 = 0
number = 1
### best hit equals best match
score = 0
best_match = "none"
best_family = "none"
method = "none"
reference = "none"
for line in f.readlines():
        line =line.strip('\t')
        token = line.split()
        token_part = token[0].split('_')
        parts = last_seq.split('_')
        if(number < int(parts[0])):
                for i in range(number,int(parts[0])):
                        print (str(i)+'\t'+'0:0:0:0'+'\t'+"none"+'\t'+"none"+'\t'+"none"+'\t'+"none")
                        x = 0
                number = int(parts[0])
        #print(line)
        ## token[1] is sequence header name
        min_len = min(int(len(PMRD_seq[token[1]])),int(token_part[2]))
        #print(min_len,int(token[3]))
        
        if( last_seq == token[0]):
                ## Case1: No mis-matches
                if( int(token[3])==int(min_len)&(int(token[4])==0)):
                        #print(int(token[3]), int(min_len),PMRD_seq[token[1]],last_seq)
                        #print("match",min_len,int(token[3]),int(PMRD_seq[token[1]]),int(token_part[2]))
                        count_0 += 1
                        if(score < 1000):
                                score = 1000
                                best_match = PMRD_seq[token[1]]
                                best_family = PMRD_family[token[1]]
                                method = PMRD_method[token[1]]
                                reference = PMRD_reference[token[1]]
                                #print(last_seq+'\t'+str(count_0)+':'+str(count_1)+':'+str(count_2)+':'+str(count_3)+'\t'+best_match+'\t'+best_family+'\t'+method+'\t'+reference)
                                #print(score)
                 ## Case2: Exactly one mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                if((int(token[3])+1-int(token[4]))==int(min_len)):
                        count_1 +=1
                        if ( score < 100):
                                best_match = PMRD_seq[token[1]]
                                best_family = PMRD_family[token[1]]
                                method = PMRD_method[token[1]]
                                reference = PMRD_reference[token[1]]
                                score = 100
                                #print(last_seq+'\t'+str(count_0)+':'+str(count_1)+':'+str(count_2)+':'+str(count_3)+'\t'+best_match+'\t'+best_family+'\t'+method+'\t'+reference)
                                #print(score)
                 ## Case3: Exactly two mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                if((int(token[3])+2-int(token[4]))==int(min_len)):
                        count_2 +=1
                        if ( score < 10):
                                best_match = PMRD_seq[token[1]]
                                best_family = PMRD_family[token[1]]
                                method = PMRD_method[token[1]]
                                reference = PMRD_reference[token[1]]
                                score = 10
                                #print(last_seq+'\t'+str(count_0)+':'+str(count_1)+':'+str(count_2)+':'+str(count_3)+'\t'+best_match+'\t'+best_family+'\t'+method+'\t'+reference)
                                #print(score)
                 ## Case4: Exactly three mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                if((int(token[3])+3-int(token[4]))==int(min_len)):
                        count_3 +=1
                        if ( score < 1):
                                best_match = PMRD_seq[token[1]]
                                best_family = PMRD_family[token[1]]
                                method = PMRD_method[token[1]]
                                reference = PMRD_reference[token[1]]
                                score = 1
                                #print(last_seq+'\t'+str(count_0)+':'+str(count_1)+':'+str(count_2)+':'+str(count_3)+'\t'+best_match+'\t'+best_family+'\t'+method+'\t'+reference)
                                #print(score)
                else:
                        continue
        else:
                print(last_seq+'\t'+str(count_0)+':'+str(count_1)+':'+str(count_2)+':'+str(count_3)+'\t'+best_match+'\t'+best_family+'\t'+method+'\t'+reference)
                count_0 = 0
                count_1 = 0
                count_2 = 0
                count_3 = 0
                best_match = "none"
                best_family = "none"
                method = "none"
                reference = "none"
                number += 1
                score = 0
        ## Case2: One mis-matches allowed ( min of both query and target seq are allowed to have one mismatch of start/end layoff
                
        
        last_seq = token[0]

"""

