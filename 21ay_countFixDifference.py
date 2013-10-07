#-----------------------------------------------------------+
#                                                           |
# 21ay_countFixDifference.py - script to calculate the ancestral alleles and fix differneces      |
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
import os,sys,getopt, re, E_get_chr_size_gff3, classGene, classMAF


### global variables
global ifile, marker, align

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
            python 21ay_countFixDifference.py -i <ifile> ### GFF3 file
                                              -m <marker> ### marker file
                                              -a <align> ### alignment file
            '''
    sys.exit(2)

### main argument to 

def options(argv):
    global ifile, marker, align
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:m:a:",["ifile=","marker=","align="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-m", "--marker"):
            marker = arg
        elif opt in ("-a", "--align"):
            align = arg
            
        
            
    logfile(ifile)
            
def hashGFF3(chrom):
    cds_coords = {}
    for line in open(ifile,'r'):
        line = line.strip()
        if len(line) > 1 and line.startswith(chrom):
            obj = classGene.GFF3(line)
            if obj.types() == 'CDS':
                for i in range(int(obj.starts()), int(obj.ends())+1, 1):
                    cds_coords[i] = ''
    return cds_coords
        
def hashAlignment(chrom, cds_coords, o2):
    objs = ''
    new_black = True
    ref_align = True
    obj_count=1
    obj = classMAF.MAF(obj_count)
    objs_list = [obj]
    for line in open(align, 'r'):
        line = line.strip()
        if line == '':
            obj_count += 1
            ref_align = True
            obj = classMAF.MAF(obj_count)
            objs_list.append(obj)
        else:
            if not line.startswith('s '):
                obj.addData(line)
            else:
                if ref_align == True:
                    obj.addRefAlign(line)
                    ref_align = False
                else:
                    obj.addtargetAlign(line)
    
    align_coords={}
    for ob in objs_list[:-1]:
        if ob.RefChro() == chrom:
            temp = ob.AlignHash()
            if ob.RefStrand() == '+':
                for i in range(len(temp)):
                    align_coords[int(ob.RefStart()) + i + 1] = temp[i+1]
                    
                    ### print mismaches
                    if temp[i+1][0] != temp[i+1][1]:
                        if (int(ob.RefStart()) + i + 1) in cds_coords:
                            o2.write(ob.RefChro()+'\t'+str(int(ob.RefStart()) + i + 1)+'\t'+temp[i+1][0]+'\t'+temp[i+1][1]+'\t'+'True'+'\n')
                        else:
                            o2.write(ob.RefChro()+'\t'+str(int(ob.RefStart()) + i + 1)+'\t'+temp[i+1][0]+'\t'+temp[i+1][1]+'\t'+'False'+'\n')
            else:
                print 'Error at line'
                print str(obj)
                sys.exit('Refernce strand is - ')
            
    return align_coords

def parseMarkers(chrom, cds_coords, align_coords, o1):
    for line in open(marker, 'r'):
        if len(line) > 0 and line.startswith(chrom):
            token = line.split('\t')
            chro = token[0]
            pos = int(token[1])
            ref = token[3]
            alt = token[4]
            string = str(chro)+'\t'+str(pos)+'\t'+str(ref)+'\t'+str(alt)
            if pos in align_coords:
                string += '\t' +  align_coords[pos][0] + '\t' + align_coords[pos][1]
            else:
                string += '\t' +  'None' + '\t' + 'None'
            
            if pos in cds_coords:
                string += '\t' + 'True'
            else:
                string += '\t' + 'False'
            
            o1.write(string+'\n')
    
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    size = E_get_chr_size_gff3.get_size(ifile)
    
    o1 = open(marker+'.medicagoAnceAllele','w')
    HEADER='Chromosome\tPosition\tRefBase\tAltBase\tLotusAlignBase\tMedicagoBase\tInCDS'
    o1.write(HEADER+'\n')
    
    o2 = open(marker+'.medicagoFixDiff','w')
    HEADER='Chromosome\tPosition\tLotusAlignBase\tMedicagoBase\tInCDS'
    o2.write(HEADER+'\n')
    
    
    
    for chrom in size:
        cds_coords = hashGFF3(chrom)
        
        ### hash all the alignment bases
        align_coords = hashAlignment(chrom, cds_coords, o2)
    
        ### pass marker list and print
        parseMarkers(chrom, cds_coords, align_coords, o1)
        
        
    o1.close()
    o2.close()
    ### close the logfile
    o.close()