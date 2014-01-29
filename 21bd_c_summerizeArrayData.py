#-----------------------------------------------------------+
#                                                           |
# 21bd_summerizeArrayData.py - script to calculate the mean and std for micro array data       |
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
import numpy as np


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
            
def parse(ifile):
    reps = [3,3,3,3,3,3,3,3]
    header = True
    for line in open(ifile, 'r'):
        line = line.strip()
        tokens = line.split('\t')
        if len(line) > 1:
            if header == True:
                string = tokens[0]
                cols = 0
                for i in range(len(reps)):
                    string += '\t' + 'SampleValues_'+ tokens[cols+1] + '\tMean_' + tokens[cols+1] +'\tStd_'+tokens[cols+1]
                    cols += reps[i]
            else:
                cols = 0
                string = tokens[0]
                for i in range(len(reps)):
                    values = [float(tokens[cols+j+1]) for j in range(reps[i])]
                    mean = round(np.mean(values),2)
                    std = round(np.std(values),2)
                    string += '\t' + '_'.join([(tokens[cols+j+1]) for j in range(reps[i])])
                    string += '\t' + str(mean) + '\t' + str(std)
                    cols += reps[i]
                    
            header = False
                            
            print string

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    
    ### parse file
    parse(ifile)
    
    ### close the logfile
    o.close()