#-----------------------------------------------------------+
#                                                           |
# 12av_removeMAFoverlaps.py - script to take longest alignment from the overlapping alignments       |
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
import os,sys,getopt, re, classMAF


### global variables
global ifile, bed
global chro
#chro = ['chr1', 'chr2' , 'chr3', 'chr4', 'chr5', 'chr6']
chro = ['chr4']
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
    global ifile, bed
    ifile = ''
    bed = False
    try:
        opts, args = getopt.getopt(argv,"hi:b",["ifile=","bed="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-b", "--bed"):
            bed = True
            
    logfile(ifile)
            
def parseMAF():
    o = open(ifile+'.nonRed','w')
    if bed == True:
        b = open(ifile+'.bed','w')
    objs = ''
    new_black = True
    ref_align = True
    obj_count=1
    obj = classMAF.MAF(obj_count)
    objs_list = [obj]
    for chrom in chro:
        for line in open(ifile, 'r'):
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
        
        chr_obj_pos = {}
        
        for ob in objs_list[:-1]:
            if ob.RefChro() == chrom:
                chr_obj_pos[int(ob.RefStart())] = ob
        
        chr_obj = []
        for key in sorted(chr_obj_pos):
            chr_obj.append(chr_obj_pos[key])
        
        
        hash_coords={}
        hash_coords_obj={}
        obj_list = {}
        obj_count = 0
        for obj in chr_obj:
            obj_count += 1
            if obj.RefChro() == chrom:
                for i in range(int(obj.RefStart()), int(obj.RefEnd())+1, 1):
                    if i in hash_coords:
                        if hash_coords_obj[i].RefSize() > obj.RefSize():
                            if obj in obj_list:
                                del obj_list[obj]
                            break
                        else:
                            hash_coords[i] = obj.RefSize()
                            if hash_coords_obj[i] in obj_list:
                                del obj_list[hash_coords_obj[i]]
                            if obj not in obj_list:
                                obj_list[obj]=''
                            hash_coords_obj[i] = obj
                    else:
                        hash_coords[i] = obj.RefSize()
                        if obj not in obj_list:
                            obj_list[obj] = ''
                        hash_coords_obj[i] = obj
                print len(obj_list)
                print chrom, '{:9,.0f}'.format(i)
                print obj.TargetChro(), obj.TargetStart()
    
        for ob in obj_list:
            o.write(ob.Data()+'\n')
            
            if bed == True:
                b.write(ob.RefChro()+'\t'+str(ob.RefStart())+'\t'+str(ob.RefEnd())+'\t'+ob.TargetChro()+'_'+str(ob.TargetStart())+'_'+str(ob.TargetEnd())+'\t'+str(ob.identity())+'\t'+str(ob.RefStrand())+'\n')
    o.close()
    if bed == True:
        b.close()
if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    
    ### parseMAF
    parseMAF()
    
    ### close the logfile
    o.close()