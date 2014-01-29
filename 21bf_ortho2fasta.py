#-----------------------------------------------------------+
#                                                           |
# 21bf_ortho2fasta.py - script to take a list of orthologs per line and a make fasta file for each ortho group     |
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
global ifile, time_start, fasta
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
            python 100b_fasta2flat.py -i <ifile> ### file with orthogroups
                                        -f <fasta> ### multi fasta file
            '''
    sys.exit(2)
    
''' orthogroup file example '''
'''
Arabidopsis	Lotus	Soybean	Medicago
AT1G11160	Lj6g2066940	PAC26331824	Medtr2g012630.1
AT2G22125	Lj6g0525130	PAC26355956	Medtr8g091470.1
AT3G18890	Lj3g3338270	PAC26292126	Medtr4g068970.1
AT1G72440	Lj3g2575990	PAC26349448	Medtr5g019130.1
'''

### main argument to 

def options(argv):
    global ifile, fasta
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:f:",["ifile=","fasta="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-f", "--fasta"):
            fasta = arg
            
    logfile(ifile)
            
def LOADfasta(file):
    first_line = True
    seq = {}
    full_headers = {}
    string = ''
    for line in open(file,'r'):
            line = line.strip()
            if len(line) > 0 :			
                    if line[0] == '>':
                            if first_line == False:
                                    if string != '': 
                                            seq[header] = string
                            string = ''
                            header = line[1:].strip().split(',')[0].split(' ')[0]
                            
                            ### special case for glycin max
                            if line[1:].startswith('Glyma'):
                                header = line[1:].split('|')[1].replace('id:','')
                            
                            full_headers[header] = line
                    else:
                            string += line
            first_line = False			
    if string != '': 
            seq[header] = string
    return seq, full_headers
                            
def print_fasta(file, seq):
    first_line = True
    for line in open(file,'r'):
        line = line.strip()
        ids = line.split()
        if first_line == False:
            out = open(ids[0]+'.fa','w')
            for id in ids:
                if (id in seq):
                    out.write('>'+id+'\n')
                    out.write(seq[id]+'\n')
                elif (id+'.1' in seq):
                    id = id+'.1'
                    out.write('>'+id+'\n')
                    out.write(seq[id]+'\n')
                elif (id+'.2' in seq):
                    id = id+'.2'
                    out.write('>'+id+'\n')
                    out.write(seq[id]+'\n')
                elif (id+'.3' in seq):
                    id = id+'.3'
                    out.write('>'+id+'\n')
                    out.write(seq[id]+'\n')  
                else:
                    print 'Error at line'
                    print line
                    sys.exit('Following id does not exit: '+str(id))
            
            out.close()
        first_line = False
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### hash the fasta file
    seq, full_headers = LOADfasta(fasta)
    
    ### print the fasta file
    print_fasta(ifile,seq)
    
    ### close the logfile
    o.close()