#-----------------------------------------------------------+
#                                                           |
# 116_runGATK.py - script to run GATK analysis              |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

# Example:
# python ~/script/python/116_runGATK.py -i 02_Stegodyphous_cdna.refined.fa.orf.tr_longest_frame


### import modules
import os,sys,getopt, re


### global variables
global bams, ref, picard, gatk, threads, variant, sort_bams

### make a logfile
import datetime
now = datetime.datetime.now()
o = open(str(now.strftime("%Y-%m-%d_%H%M."))+'logfile','w')



### write logfile

def logfile(infile):
    o.write("Program used: \t\t%s" % "116_runGATK.py"+'\n')
    o.write("Program was run at: \t%s" % str(now.strftime("%Y-%m-%d_%H%M"))+'\n')
    o.write("Infile used: \t\t%s" % infile+'\n')
            
    
def help():
    print '''
            python 116_runGATK.py
                                -b <bams> [One bam file per sample seperated by commas]
                                -r <ref> [Reference sequence]
                                -p <picard> [Path to picard folder MarkDuplicates.jar]
                                -g <gatk> [Path to GATK folder containing GenomeAnalysisTK.jar]
                                -t <threads> [Number of threads to be used]
                                
            
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    
    global bams, ref, picard, gatk, threads, variant, sort_bams
    
    bams = ''
    ref = ''
    picard = ''
    gatk = ''
    threads = 1
    variant = ''
    sort_bams = False
    
    try:
        opts, args = getopt.getopt(argv,"hb:r:p:g:t:v:s",["bams=","ref=","picard=","gatk=","threads=","variant=", "sort_bams="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-b", "--bams"):
            bams = arg
        elif opt in ("-r", "--ref"):
            ref = arg
        elif opt in ("-p", "--picard"):
            picard = arg
        elif opt in ("-g", "--gatk"):
            gatk = arg
        elif opt in ("-t", "--threads"):
            threads = str(arg)
        elif opt in ("-v", "--variant"):
            variant = arg
        elif opt in ("-s", "--sort"):
            sort_bams = True
            
            
    logfile(bams)

def Index():
    if not os.path.isfile(ref + '.fai'):
        os.system('nice -n 19 samtools faidx '+ ref)
                              
                          
def files():
    files = []
    print bams
    for file in bams.split(','):
        print file
        files.append(file.strip())
    return files

def sortBams(file_list):
    if sort_bams == True:
        file_list_bams = []
        for file in file_list:
            os.system('samtools view -bt ' + ref +'.fai ' + file +' | samtools sort - ' + file + '_sorted')
            file_list_bams.append(file+'_sorted.bam')
            
        return file_list_bams
    return file_list
    
def MarkDuplicates(file_list):
    for file in file_list:
        print 'Marking duplicates for', file
        os.system('java -jar '+picard+'/MarkDuplicates.jar ' +' INPUT='+ file + ' OUTPUT=' + file+'.dedup.bam' +' METRICS_FILE='+ file +'.dups')
        
    
def ReAlign(file_list):
    for file in file_list:
        print 'Realigning', file
        os.system('java -jar '+gatk+'/GenomeAnalysisTK.jar -U -T RealignerTargetCreator '+' -I ' + file+'.dedup.bam' +' -nt '+ threads +' -R ' + ref +' -o '+ file+'.intervals')
        os.system('java -jar '+gatk+'/GenomeAnalysisTK.jar -U -T IndelRealigner ' + ' -targetIntervals ' + file+'.intervals '+ ' -I ' + file+'.dedup.bam' +' -R ' + ref \
                  + ' -o ' + file+'.realigned.bam')
    
        
def UnifiedGenotyper(file_list):
    if variant == '':
        in_string = ''
        ### make input string
        for file in file_list:
            in_string += ' -I '+file
        print 'Running UnifiedGenotyper'
        os.system(' java -jar '+gatk+'/GenomeAnalysisTK.jar '\
        + ' -R ' + ref \
        + ' -T  UnifiedGenotyper '\
        + in_string \
        + ' -nt ' + threads \
        + ' -o snps.90.raw.vcf ' \
        + '-stand_call_conf 90 ' \
        + '-stand_emit_conf 10.0 '\
        + '-dcov 200 ')
        
def recal(file_list):
    global variant
    if variant == '':
        variant = 'snps.raw.90.vcf'
    for file in file_list:
        print 'Running BaseRecalibrator for ', file
        os.system('java -jar '+gatk+'/GenomeAnalysisTK.jar -U -T BaseRecalibrator ' + ' -knownSites ' +variant + ' -I ' + file+'.realigned.bam '+ ' -R ' + ref \
              + ' -o ' + file+'.recal.bam')

def ReduceReads(file_list):
    for file in file_list:
        print 'Running ReduceReads for ', file
        os.system('java -jar '+gatk+'/GenomeAnalysisTK.jar -U -T ReduceReads -I ' + file+'.realigned.bam '+ ' -R ' + ref \
              + ' -o ' + file+'.reduced.bam')
        
def ReUnifiedGenotyper(file_list):
    print 'Running ReUnifiedGenotyper'
    if variant == '':
        in_string = ''
        ### make input string
        for file in file_list:
            in_string += ' -I '+file
        
        os.system(' java -jar '+gatk+'/GenomeAnalysisTK.jar '\
        + ' -R ' + ref \
        + ' -T  UnifiedGenotyper '\
        + in_string \
        + ' -nt ' + threads \
        + ' -o snps.raw.vcf ')

def BuildErrorModelWithVQSR(file , var):
    os.system('java -jar '+gatk+'/GenomeAnalysisTK.jar '\
    + ' -T VariantRecalibrator ' \
    + ' -R '+ ref \
    + ' -input '+ file \
    + ' -recalFile output.recal ' \
    + ' -tranchesFile output.tranches ' \
    + ' -nt ' + threads \
    + ' -mode ' + var)

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    ### check if index exits
    Index()
    
    
    ### return the list of the bam files
    file_list = files()

    ### sort Bams file
    file_list = sortBams(file_list)
    
    ### mark duplicates
    MarkDuplicates(file_list)
    
    ### realign the reads
    ReAlign(file_list)
    
    ### call UnifiedGenotyper to make a primary list of variants
    UnifiedGenotyper(file_list)
    
    ### Baserecalibration
    recal(file_list)
    
    ### reducing BAM files
    ReduceReads(file_list)
    
    ### Run UnifiedGenotyper
    ReUnifiedGenotyper(file_list)
    
    
    ### BuildErrorModelWithVQSR
    #BuildErrorModelWithVQSR('snps.raw.vcf', 'SNP')
    #BuildErrorModelWithVQSR('snps.raw.vcf', 'INDEL')
    
    ### close the logfile
    o.close()