'''
Created on May 03, 2013

@author: vgupta
'''

### 28c_gff3_validator.py - /Users/vikas0633/Desktop/script/python/ - Script to validate a gff3 file

### this script is made to test the RAP-DB GFF file set
### A file must meet the following pre-requirements
'''
1. File has at least two types of the features, mRNA and CDS
2. File has be in the blocks where CDS follow corresponding mRNA
3. each mRNA must have a locus ID
4. CDS must have parent ID
'''
### Sample
'''
chr01    irgsp1_prediction    mRNA    69675    70131    -1.0735    +    .    ID=Os01t0101500-00;Name=Os01t0101500-00;Locus_id=Os01g0101500
chr01    irgsp1_prediction    CDS    69675    69962    -1.96261    +    0    Parent=Os01t0101500-00
chr01    irgsp1_prediction    CDS    70027    70131    -1.97213    +    0    Parent=Os01t0101500-00
chr01    irgsp1_prediction    mRNA    220411    220776    -1.98506    -    .    ID=Os01t0104300-00;Name=Os01t0104300-00;Locus_id=Os01g0104300
chr01    irgsp1_prediction    CDS    220411    220776    -1.98506    -    0    Parent=Os01t0104300-00
chr01    irgsp1_prediction    mRNA    289007    289294    -0.294654    -    .    ID=Os01t0105450-00;Name=Os01t0105450-00;Locus_id=Os01g0105450
chr01    irgsp1_prediction    CDS    289007    289294    -0.294654    -    0    Parent=Os01t0105450-00 
'''

### Script will test test for following conditions
'''
1. Each gene model must have only one Locus_id
2. Gene model ID must correspond to Lotus_id
3. Any non-mRNA feature type must have corresponding parent
4. Any non-mRNA feature type must not have co-ordinates laying outside the mRNA boundries
5. For any feature type start must be smaller than or equal to end
6. A feature type should not have an overlapping block with another block of same feature   
'''

import os,sys,getopt, re



### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    

### main argument to 

def options(argv):
    infile = ''
    gff3 = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print '''
            python 28c_gff3_validator.py -i <inputfile>
            '''
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print '''
                python 28c_gff3_validator.py -i <inputfile>
                '''
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
            
    logfile(infile)
            
    return infile


def process_mRNA(line, token):
    
    ID = line.split('ID=')[1].split(';')[0]
    Locus_id = line.split('Locus_id=')[1].split(';')[0]
    
    ### check for condition 1 that there is only one Locus ID
    if len(line.split('Locus_id=')) >2:
        print sys.exit("Error: Multiple Locus_id associated with the mRNA " + "on line\n" + line ) 

    ### check for condition 2 that ID corresponds to Locus_id
    if Locus_id != ID.split('-')[0].replace('t','g'):
        print sys.exit("Error: Unexpected ID or Locus_id values on line\n"+ line)
    
    return int(token[3]), int(token[4]), ID


def count_features(feature_counts,feature):
    if feature not in feature_counts:
        feature_counts[feature] = 1
    else:
        feature_counts[feature] = 1
    
    return feature_counts
                
def process_feature(line, token, features, ID):
    feature = token[2]

    parent_id = line.split("Parent=")[1].split(';')[0]
    
    ### check if the feature has the parent
    if len(line.split("Parent=")) != 2:
        print sys.exit("Error: "+ feature+ " has unexpected number of parents at line\n"+ line )
        
    ### check condition 3, if feature has correct parent
    if parent_id != ID:
        print sys.exit("Error: "+ feature + " has wrong parent on the line\n" + line)
        
    if feature not in features:
        features[feature] = {}
        for i in range(int(token[3]),int(token[4])+1):
                    features[feature][i]=''
    else:
        coords = features[feature]
        
        ### check condition 6 that blocks doesn't overlap 
        if (int(token[3]) in coords) or (int(token[4]) in coords):
            print sys.exit("Error: "+ feature + " feature type has overlapping block at line\n"+line)
        
        ### add the coord to the block
        for i in range(int(token[3]),int(token[4])+1):
                    features[feature][i]=''
        
    return features

### function to check condition 5 where start must be before end
def test_coords(token,line):
    
    if int(token[3]) > int(token[4]):
        print sys.exit("Error: start greater than end at line\n"+line)

### function to test condition 4 where no feature extends over mRNA boundries
def test_boundries(features, st, en, ID):
    for feature in features:
        
        if (min(features[feature]) < st) or (max(features[feature]) > en):
            print sys.exit("Error: "+feature+" extends beyond mRNA co-ordinate range at mRNA ID "+ID)
    
### process the GFF3 file block by block
def parse_file(file):
    
    ### define variables for the block 
    new_block = False
    feature_counts = {}
    for line in open(file,'r'):
        line = line.strip()
        token = line.split('\t')
        feature = token[2]
        
        test_coords(token,line)
        
        if len(line) > 1:
            if line[0]!='!' and line[0]!='#':
                
                ### counts the features
                feature_counts = count_features(feature_counts, feature)
                
                if feature == "mRNA":
                    if new_block == True:
                        test_boundries(features, mRNA_st, mRNA_en, ID)
                    features = {}
                    mRNA_st, mRNA_en, ID = process_mRNA(line,token)
                    new_block = True
                
                else: 
                    features = process_feature(line, token, features, ID)   
                
         


if __name__ == "__main__":
    
    file = options(sys.argv[1:])
    
    ### parse the GFF3 file
    parse_file(file)
    
    ### close the logfile
    o.close()