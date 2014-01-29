#-----------------------------------------------------------+
#                                                           |
# 21be_MakeHeatMap.py - script to make heatMaps with array data |
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
global ifile, time_start
global GeneList, SampleList, ArrayData, BlastResults, UseAvg
start_time = time.time()

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')

### return time eclapsed
def PrinteclapsedTime():
    diff = time.time() - start_time
    minutes, seconds = int(diff)/60, diff % 60
    print('Time taken Min:Sec ==> ' + str(minutes) + ':' + str(round(seconds,2))+'\n'+'\n')


### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "21be_MakeHeatMap.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 21be_MakeHeatMap.py -g <GeneList> ### one gene per row
                                        -s <SampleList> ### column names in the Array data table
                                        -d <ArrayData> ### table with all the array data
                                        -b <BlastResults> ### correspondance table between GeneID and ProbeID
                                        -i <IndProbe> ### use this option to plot individual probes
                                        -a <anno> ### tab delimited file geneID and annotations as two colums
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global GeneList, SampleList, ArrayData, BlastResults, UseAvg, anno
    UseAvg = True
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hg:s:d:b:ia:",["GeneList=","SampleList=","ArrayData=","BlastResults=","IndProbe","anno="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-g", "--GeneList"):
            GeneList = arg
        elif opt in ("-s", "--SampleList"):
            SampleList = arg
        elif opt in ("-d", "--ArrayData"):
            ArrayData = arg
        elif opt in ("-b", "--BlastResults"):
            BlastResults = arg
        elif opt in ("-i", "--IndProbe"):
            UseAvg  = False
        elif opt in ("-a", "--anno"):
            anno  = arg
            
    logfile(ifile)
            

def hashArray():
    ArrayHash = {}
    header = True
    for line in open(ArrayData, 'r'):
        line = line.strip()
        tokens = line.split('\t')
        if len(line) > 0 and not line.startswith('#'): 
            if header == True:
                samples = tokens[1:]
            else:
                ArrayHash[tokens[0]] = tokens[1:]
        header = False
    
    return samples, ArrayHash

def hashBlast():
    BlastHash = {}
    header = True
    for line in open(BlastResults, 'r'):
        line = line.strip()
        tokens = line.split('\t')
        if len(line) > 0 and not line.startswith('#'): 
            if header == True:
                samples = tokens[1:]
            else:
                if tokens[2] in BlastHash:
                    BlastHash[tokens[2]].append(tokens[0])
                else:
                    BlastHash[tokens[2]] = [tokens[0]]
            header = False
    
    return BlastHash

def hash_anno():
    AnnoHash = {}
    for line in open(anno, 'r'):
        line = line.strip()
        tokens = line.split('\t')
        g_id = '.'.join(tokens[0].split('.')[:-1])
        AnnoHash[g_id] = tokens[1].replace(' ','_')
    
    return AnnoHash

def printOut(samples, ArrayHash, BlastHash, AnnoHash):
    SampleDict = {}
    '''
    for line in open(SampleList, 'r'):
        line = line.strip()
        SampleDict[line.split()[0]] = ''
    '''
    ### print output to a file
    o = open(GeneList+'.arrayData', 'w')
    header = '\t'.join([samples[3*i+1] for i in range(len(samples)/3)])
    o.write('ProbeID\t'+'Annotation\t'+header+'\n')
    for line in open(GeneList, 'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            gene_id = line.split()[0]
            
            if gene_id in BlastHash:
                if UseAvg == False:
                    for i in range(len(BlastHash[gene_id])):
                        probe_id = BlastHash[gene_id][i]
                        #o.write(gene_id+'.'+str(i+1)+'_'+probe_id+'\n')
                        o.write(gene_id+'.'+str(i+1)+'_'+probe_id+'\t'+'\t'.join(ArrayHash[probe_id])+'\n')
                
                else:
                    high_std = False
                    probe_id = BlastHash[gene_id][0]
                    avg = [0 for i in range(len(ArrayHash[probe_id])/3)]
                    for i in range(len(BlastHash[gene_id])):
                        probe_id = BlastHash[gene_id][i]
                        for j in range(len(avg)):
                            
                            ### check for standard deviation
                            #if float(ArrayHash[probe_id][3*j+2])/float(ArrayHash[probe_id][3*j+1]) < 0.2:
                            avg[j] += float(ArrayHash[probe_id][3*j+1])
                            
                        #print ArrayHash[probe_id][3*0+1], avg[0], len(BlastHash[gene_id])
                    avg = [str(round(avg[i]/float(len(BlastHash[gene_id])),2)) for i in range(len(avg))]
                    o.write(gene_id+'\t'+AnnoHash[gene_id]+'\t'+'\t'.join(avg)+'\n')
            else:
                print "GeneID not found: ", gene_id 
    o.close()
    
    
if __name__ == "__main__":
    
    options(sys.argv[1:])
    print 'Hold on reading data: '
    PrinteclapsedTime()
    
    ### hash the array data
    samples, ArrayHash = hashArray()
    print 'Finished loading array data: '
    PrinteclapsedTime()
    
    ### hash the blast results - correspondance table
    BlastHash = hashBlast()
    print 'Finished loading blast data: '
    PrinteclapsedTime()
    
    ### hash annotations
    AnnoHash = hash_anno()
    
    ### print the gene based expressions
    printOut(samples, ArrayHash, BlastHash, AnnoHash)
    
    print 'Finished printing gene-wise table: '
    PrinteclapsedTime()
    
    
    ### run the Rscript
    #os.system('R --vanilla < ~/script/R/14_heatMap.R')
    
    ### close the logfile
    o.close()