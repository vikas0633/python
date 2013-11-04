#-----------------------------------------------------------+
#                                                           |
# 119_vcfParser.py - script to parse vcf format file       |
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
import os,sys,getopt, re, classVCF, time


### global variables
global ifile

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
            

def calc_dist(g1, g2):
    dist = 0
    
    if g1 == '0/0' and g2 == '0/1':
        dist += 0.5
    elif g1 == '0/1' and g2 == '0/0':
        dist += 0.5 
    elif g1 == '0/1' and g2 == '1/1':
        dist += 0.5
    elif g1 == '1/1' and g2 == '0/1':
        dist += 0.5
    
    elif g1 == '0/0' and g2 == '1/1':
        dist += 1
    elif g1 == '1/1' and g2 == '0/0':
        dist += 1
    
    return dist
        
        

def printOut(dist_mat, genotypes):
    o = open(ifile+'.dist','w')
    o.write(str(len(genotypes)))
    for i in range(len(genotypes)):
        o.write('\n'+genotypes[i])
        for j in range(len(genotypes)):
            o.write('\t'+str(dist_mat[i,j]))
    o.close()

def parse():
    count = 0
    then = time.time()
    for line in open(ifile, 'r'):
        if len(line) > 0 and not line.startswith('##'):
            line = line.strip()
            obj = classVCF.VCF(line)
            
            
            count += 1
            if count%10000 == 0:
                diff = time.time() - then
                minutes, seconds = int(diff)/60, diff % 60
                print 'Number of markers processed: ', '{:9,.0f}'.format(count)
                print('Time taken Min:Sec ==> ' + str(minutes) + ':' + str(round(seconds,2)))
            
            if line.startswith('#'):
                genotypes = obj.genotypes()
                g_count = len(genotypes)
                dist_mat = {}
                for i in range(g_count):
                    for j in range(g_count):
                        dist_mat[i,j] = 0
            else:
                for i in range(g_count):
                    for j in range(g_count):
                        geno1 = obj.genotype(i)
                        geno2 = obj.genotype(j)
                        dist_mat[i,j] += calc_dist(geno1, geno2)
    
    printOut(dist_mat, genotypes)

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### parse vcf
    parse()
    
    
    ### close the logfile
    o.close()