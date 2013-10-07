#-----------------------------------------------------------+
#                                                           |
# 21ba_getGeneBasedAlign.py - script to extract gene specific information from MAF file      |
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
global ifile, candidate, align

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
    global ifile, candidate, align
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:a:c:",["ifile=","candidate="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-a", "--align"):
            align = arg
        elif opt in ("-c", "--candidate"):
            candidate = arg
        
            
        
            
    logfile(ifile)

def hashAlignment(chrom, cds_coords, start):
    objs = ''
    new_black = True
    ref_align = True
    print_flag = False
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
    
    query = candidate+'_Lj\t'
    target = candidate+'_Mt\t'
    blocks = {}
    align_coords={}
    for ob in objs_list[:-1]:
        if ob.RefChro() == chrom:
            temp = ob.AlignHash()
            if ob.RefStrand() == '+':
                for i in range(len(temp)):
                    if (int(ob.RefStart()) + i) in cds_coords:
                        if len(blocks) > 1 and last_block_len != len(blocks):
                            query += str('-'*99)
                            target += str('-'*99)
                            print_flag = True
                        query += temp[i][0]
                        target += temp[i][1]
                        last_block_len = len(blocks)
                        blocks[ob] = ''
                        
                        
            else:
                print 'Error at line'
                print str(obj)
                sys.exit('Refernce strand is - ')
    if print_flag == True:
        print query
        print target
    
def hashGFF3(chrom):
    cds_coords = {}
    cds_bound = {}
    for line in open(ifile,'r'):
        line = line.strip()
        if len(line) > 1 and line.startswith(chrom):
            obj = classGene.GFF3(line)
            if obj.types() == 'CDS' and obj.get_parent() == candidate:
                cds_bound[int(obj.starts())] = ''
                cds_bound[int(obj.ends())] = ''
                
                for i in range(int(obj.starts()), int(obj.ends())+1, 1):
                    cds_coords[i] = ''
                    

    
    if len(cds_coords) != 0:
        if obj.strands() == '+':
            start = min(cds_bound)
        else:
            start = max(cds_bound)
        hashAlignment(chrom, cds_coords, start)
        
    
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    size = E_get_chr_size_gff3.get_size(ifile)
    
    
    
    for chrom in size:
        cds_coords = hashGFF3(chrom)
        
    
        
        
    ### close the logfile
    o.close()