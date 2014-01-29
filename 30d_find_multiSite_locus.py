#-----------------------------------------------------------+
#                                                           |
# 30d_find_multiSite_locus.py - script to find the locus with multiple significant degradome sites      |
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
import os,sys,getopt, re, time, numpy
chro_list = ['chr0','chr1','chr2','chr3','chr4','chr5','chr6','chloro','mito']

### global variables
global ifile, time_start
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
    o.write("Program used: \t\t%s" % "30d_find_multiSite_locus.py "+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 30d_find_multiSite_locus.py -i <ifile>
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
            
    logfile(ifile)
            

def find_sites(chro):
    new_block = True
    last_pos = 0
    SignList = []
    for line in open(ifile, 'r'):
        line = line.strip()
        tokens = line.split('\t')
        pos = int(tokens[1])
        sign = float(tokens[4])
        if tokens[0] == chro:
            if pos - last_pos > 1000 :
                if numpy.sum(SignList) > 4 and len(SignList) >= 4 and numpy.mean(SignList)>0.2:
                    print chro, last_pos
                avg_sign = 0
                SignList = []
                new_block = True
            else:
                SignList.append(sign)
                new_block = False
        
            last_pos = pos
        

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### run one chromosome at a time
    for chro in chro_list:
        print 'Processing the :', chro
        PrinteclapsedTime()
        
        find_sites(chro)
    

    ### find the sites 
    
    
    
    
    ### close the logfile
    o.close()