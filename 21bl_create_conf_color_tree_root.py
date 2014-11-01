#-----------------------------------------------------------+
#                                                           |
# template.py - Python template                             |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                        |
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
import os,sys,getopt, re

import threading

from multiprocessing import Process, Queue, Manager
from threading import Thread
import classGene
### global variables
global infile

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')


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
    size_sorted={}
    for w in sorted(size, key=size.get, reverse=False):
        size_sorted[w]=size[w]
    return size_sorted


### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "100b_fasta2flat.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 100b_fasta2flat.py -i <infile>
            '''
    sys.exit(2)

def temp(file):
    return

### main argument to 

def options(argv):
    global infile, threads
    infile = ''
    threads = 2
    
    try:
        opts, args = getopt.getopt(argv,"hi:t:",["infile=","threads="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-t", "--threads"):
            threads = int(arg)
            
    


def convert_to_rgb(minval, maxval, val, colors):
    max_index = len(colors)-1
    v = float(val-minval) / float(maxval-minval) * max_index
    i1, i2 = int(v), min(int(v)+1, max_index)
    (r1, g1, b1), (r2, g2, b2) = colors[i1], colors[i2]
    f = v - i1
    return int(r1 + f*(r2-r1)), int(g1 + f*(g2-g1)), int(b1 + f*(b2-b1))

def max_min():
    global min_value, max_value
    max_value = -1000000000
    min_value = 1000000000
    for line in open(infile,'r'):
        line = line.strip()
        tokens = line.split('\t')
        if tokens[2] != 'NULL':
            value = float(tokens[2])
            if min_value > value:
                min_value = value
            elif value > max_value:
                max_value = value
    
    return min_value, max_value

def print_rgb():
    at_up = 0
    at_down = 0
    lj_up = 0
    lj_down = 0
    gm_up = 0
    gm_down = 0
    mt_up = 0
    mt_down = 0
    os_up = 0
    os_down = 0
    zm_up = 0
    zm_down = 0
    la_up = 0
    la_down = 0
    br_up = 0
    br_down = 0
    colors = [(0, 0, 255), (0, 255, 0), (255, 0, 0)]
    col = [(0,0,255),(192,192,192),(255,0,0),(0,128,255),(255,102,102)]
    for line in open(infile,'r'):
        line = line.strip()
        tokens = line.split('\t')
        if tokens[2] != 'NULL':
            value = float(tokens[2])
            
            ### Arabidopsis
            if tokens[1].startswith('Arabidopsis'):
                if value <= 3 :
                    r, g, b  = col[0]
                    at_down += 1
                elif 3 <= value < 7:
                    r, g, b  = col[1]
                elif value >= 7:
                    r, g, b  = col[2]
                    at_up += 1
                print 'prefix\t' + tokens[0] + '\t' + str(r)+' '+str(g)+' '+str(b)+'\t'+'\t'+ '2'
            ### Lotus
            if tokens[1].startswith('Lotus'):
                if value <= 3 :
                    r, g, b  = col[0]
                    lj_down += 1
                elif 3 <= value < 4:
                    r, g, b  = col[3]
                elif 4 <= value < 5.5:
                    r, g, b  = col[1]
                elif 5.5 <= value < 6:
                    r, g, b  = col[4]
                elif value >= 6:
                    r, g, b  = col[2]
                    lj_up += 1
                print 'prefix\t' + tokens[0] + '\t' + str(r)+' '+str(g)+' '+str(b)+'\t'+'\t'+ '2'
            ### Soybean
            if tokens[1].startswith('Soybean'):
                if value <= 1 :
                    r, g, b  = col[0]
                    gm_down += 1
                elif 1 <= value < 1.5:
                    r, g, b  = col[3]
                elif 1.5 <= value < 2.0:
                    r, g, b  = col[1]
                elif 2.0 <= value < 2.5:
                    r, g, b  = col[4]
                elif value >= 2.5:
                    r, g, b  = col[2]
                    gm_up += 1
                print 'prefix\t' + tokens[0][:-1] + '\t' + str(r)+' '+str(g)+' '+str(b)+'\t'+'\t'+ '2'
            ### Medicago
            if tokens[1].startswith('Medicago'):
                if value <= 0.5 :
                    r, g, b  = col[0]
                    mt_down += 1
                elif 0.5 <= value < 1.5:
                    r, g, b  = col[3]
                elif 1.5 <= value < 1.75:
                    r, g, b  = col[1]
                elif 1.75 <= value < 2:
                    r, g, b  = col[4]
                elif value >= 2:
                    r, g, b  = col[2]
                    mt_up += 1
                print 'contain\t' + tokens[0].split('|')[1][:-1] + '\t' + str(r)+' '+str(g)+' '+str(b)+'\t'+'\t'+ '2'
            ### Rice
            if tokens[1].startswith('Rice'):
                if value <= 2.5 :
                    r, g, b  = col[0]
                    os_down += 1
                elif 2.5 <= value < 6.0:
                    r, g, b  = col[3]
                elif 6.0 <= value < 7.5:
                    r, g, b  = col[1]
                elif 7.5 <= value < 9:
                    r, g, b  = col[4]
                elif value >= 9:
                    r, g, b  = col[2]
                    os_up += 1
                print 'prefix\t' + tokens[0][:-1] + '\t' + str(r)+' '+str(g)+' '+str(b)+'\t'+'\t'+ '2'
            ### Maize
            if tokens[1].startswith('Maize'):
                if value <= 6 :
                    zm_down += 1
                    r, g, b  = col[0]
                elif 6 <= value < 8:
                    r, g, b  = col[3]
                elif 8 <= value < 10:
                    r, g, b  = col[1]
                elif 10 <= value < 13:
                    r, g, b  = col[4]
                elif value >= 13:
                    r, g, b  = col[2]
                    zm_up += 1
            
                #r, g, b = convert_to_rgb(min_value, max_value, value, colors)
                print 'prefix\t' + tokens[0] + '\t' + str(r)+' '+str(g)+' '+str(b)+'\t'+'\t'+ '2'
                
            ### Lupin
            if tokens[1].startswith('Lupin'):
                if value <= -1 :
                    la_down += 1
                    r, g, b  = col[0]
                elif -1 <= value < 0:
                    r, g, b  = col[3]
                elif 0 <= value < 1.5:
                    r, g, b  = col[1]
                elif 1.5 <= value < 2:
                    r, g, b  = col[4]
                elif value >= 2:
                    r, g, b  = col[2]
                    la_up += 1
            
                #r, g, b = convert_to_rgb(min_value, max_value, value, colors)
                print 'prefix\t' + tokens[0] + '\t' + str(r)+' '+str(g)+' '+str(b)+'\t'+'\t'+ '2'
            
            ### Brassica
            if tokens[1].startswith('Brassica'):
                if value <= -5 :
                    la_down += 1
                    r, g, b  = col[0]
                elif -5 <= value < -2:
                    r, g, b  = col[3]
                elif -2 <= value < 1:
                    r, g, b  = col[1]
                elif 1 <= value < 5:
                    r, g, b  = col[4]
                elif value >= 5:
                    r, g, b  = col[2]
                    la_up += 1
            
                #r, g, b = convert_to_rgb(min_value, max_value, value, colors)
                print 'prefix\t' + tokens[0] + '\t' + str(r)+' '+str(g)+' '+str(b)+'\t'+'\t'+ '2'
                

    print >> sys.stderr, "at_down :"+str(at_down)
    print >> sys.stderr, "at_up :"+str(at_up)
    print >> sys.stderr, "lj_down :"+str(lj_down)
    print >> sys.stderr, "lj_up :"+str(lj_up)
    print >> sys.stderr, "mt_down :"+str(mt_down)
    print >> sys.stderr, "mt_up :"+str(mt_up)
    print >> sys.stderr, "gm_down :"+str(gm_down)
    print >> sys.stderr, "gm_up :"+str(gm_up)
    print >> sys.stderr, "os_down :"+str(os_down)
    print >> sys.stderr, "os_up :"+str(os_up)
    print >> sys.stderr, "zm_down :"+str(zm_down)
    print >> sys.stderr, "zm_up :"+str(zm_up)

if __name__ == "__main__":
    

    options(sys.argv[1:])
        
    start_time = datetime.datetime.now()
    print >> sys.stderr, "Running temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Input count: " + str(temp(infile))
    
    max_min()
    
    print_rgb()
    
    print >> sys.stderr, "Output count: " + str(temp(infile))
    print >> sys.stderr, "Completed temp script: " + str(datetime.datetime.now())
    print >> sys.stderr, "Time take to complete: " + str(datetime.datetime.now() - start_time)
