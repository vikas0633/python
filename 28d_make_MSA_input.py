#-----------------------------------------------------------+
#                                                           |
# 28d_make_MSA_input.py - script to make input for Multiple Sequence alignment      |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
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
import os,sys,getopt, re, time


### global variables
global ifile, time_start, fasta, cyc
start_time = time.time()

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')

### return time eclapsed
def PrinteclapsedTime():
    diff = time.time() - start_time
    minutes, seconds = int(diff)/60, diff % 60
    print('Time taken Min:Sec ==> ' + str(minutes) + ':' + str(round(seconds,2)))


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
    global ifile, fasta, cyc
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:f:p:",["ifile=","fasta=", "cyc="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-f", "--fasta"):
            fasta = arg
        elif opt in ("-p", "--cyc"):
            cyc = arg
    logfile(ifile)
            

def hashVCFasta():
    hash = {}
    for line in open(fasta,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            if line.startswith('>'):
                header = line[1:]
                hash[header] = ''
            else:
                hash[header] += line
    return hash
                            
def hashPlantFasta():
    hash = {}
    hash_full_header = {}
    for line in open(cyc,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            if line.startswith('>'):
                header = line.split('|')[2].strip()
                full_header = line[1:]
                if header in hash_full_header:
                    hash_full_header[header].append(full_header)
                else:
                    hash_full_header[header] = [full_header]
                hash[full_header] = ''
            else:
                hash[full_header] += line
    return hash, hash_full_header

def hashCorr():
    hash = {}
    for line in open(ifile,'r'):
        line = line.strip()
        tokens = line.split('\t')
        if tokens[7] not in hash:
            hash[tokens[7]] = [tokens[1]]
        else:
            hash[tokens[7]].append(tokens[1])
    return hash

def printOut(hash_vc_fasta, hash_plant_fasta, hash_corr, hash_full_header):
    count = 0
    for key in sorted(hash_corr):
        count += 1
        o = open(str(count)+'.mfa','w')
        try:
            for t_id in hash_corr[key]:
                o.write('>'+t_id+'\n')
                o.write(hash_vc_fasta[t_id]+'\n')
            for p_id in hash_full_header[key]:
                o.write('>' + p_id + '\n')
                o.write(hash_plant_fasta[p_id]+'\n')
        except:
            continue
        o.close()

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### hash vc fasta file
    hash_vc_fasta = hashVCFasta()
    
    ### hash the plantcyc fasta file
    hash_plant_fasta, hash_full_header = hashPlantFasta()
    
    ### hash the correspondance file
    hash_corr = hashCorr()
    
    ### write the output
    printOut(hash_vc_fasta, hash_plant_fasta, hash_corr, hash_full_header)
    
    ### close the logfile
    o.close()