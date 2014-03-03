#-----------------------------------------------------------+
#                                                           |
# 31a_reformat_gff3.py - script to replace the ref column of   |
#                       gff3 by priority                    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 19/02/2014                                       |
# UPDATED: 19/02/2014                                       |
#                                                           |
# DESCRIPTION:                                              | 
# script to replace the ref column of gff3 by priority      |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# nice -n 19 python -i $file.gff3 -p priority.csv > `pwd`"/gene_evidences/""${file_name[@]:(-1)}"


### import modules
import os,sys,getopt, re
import classGene

### global variables
global infile, prior

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



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

### main argument to 

def options(argv):
    global infile, prior
    infile = ''
    prior = ''
    try:
        opts, args = getopt.getopt(argv,"hi:p:",["ifile=","priority="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-p", "--priority"):
            prior = arg
            
    logfile(infile)
    return infile
    

def hash_prior():
    first_line = True
    list_file_names = []
    priority = []
    for line in open(prior, 'r'):
        line = line.strip()
        if first_line == True:
            names = line.split(',')
            for name in names:
                list_file_names.append(name)
            index = list_file_names.index((infile.split('/')[-1]).replace('.gff3',''))
            first_line = False
        else:
            prior_n = index + 1
    return prior_n

def change_gff3(prior_n):
    for line in open(infile, 'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            obj = classGene.GFF3(line)
            
            string = obj.seqids() + '\t' + \
            str(prior_n) + '\t' + \
            obj.types() + '\t' + \
            obj.starts() + '\t' + \
            obj.ends() + '\t' + \
            obj.scores() + '\t' + \
            obj.strands() + '\t' + \
            obj.phases() + '\t' + \
            obj.attributes()
        
            print string 

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### hash the priority file
    prior_n = hash_prior()
    
    ### change the gff3 file
    change_gff3(prior_n)
    
    ### close the logfile
    o.close()