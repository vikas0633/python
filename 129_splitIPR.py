#-----------------------------------------------------------+
#                                                           |
# 129_splitIPR.py - Script to split IPRscan output          |
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
import os,sys,getopt, re


### global variables
global ifile, out

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
    global ifile, out
    ifile = ''
    try:
        opts, args = getopt.getopt(argv,"hi:o:",["ifile=","out="])
    except getopt.GetoptError:
        help()
    for opt, arg in opts:
        if opt == '-h':
            help()
        elif opt in ("-i", "--ifile"):
            ifile = arg
        elif opt in ("-o", "--out"):
            out = arg
            
    logfile(ifile)
    
def IPRstring():
    string1 = '<EBIInterProScanResults   > \n' + \
'<Header> \n' + \
'	<program name="InterProScan" version="4.8" citation="PMID:11590104"/> \n' + \
'	<parameters> \n' + \
'		<sequences total="1"/> \n' + \
'		<databases total="13"> \n' + \
'			<database number="1" name="PRODOM" type="sequences"/> \n' + \
'			<database number="2" name="PRINTS" type="matrix"/> \n' + \
'			<database number="3" name="PIR" type="model"/> \n' + \
'			<database number="4" name="PFAM" type="model"/> \n' + \
'			<database number="5" name="SMART" type="model"/> \n' + \
'			<database number="6" name="TIGRFAMs" type="model"/> \n' + \
'			<database number="7" name="PROFILE" type="strings"/> \n' + \
'			<database number="8" name="PROSITE" type="strings"/> \n' + \
'			<database number="9" name="SUPERFAMILY" type="model"/> \n' + \
'			<database number="10" name="GENE3D" type="model"/> \n' + \
'			<database number="11" name="PANTHER" type="model"/> \n' + \
'			<database number="12" name="SIGNALP" type="model"/> \n' + \
'			<database number="13" name="TMHMM" type="model"/> \n' + \
'		</databases> \n' + \
'	</parameters> \n' + \
'</Header>\n'
    
    string2 = '</EBIInterProScanResults> \n'
    
    return string1, string2
            
def readIn(string1, string2):
    print_flag = False
    i = 0
    for line in open(ifile,'r'):
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            if line.startswith('<protein id='):
                i += 1
                g_id = line.split('"')[1]
                if i%100 == 0:
                    print g_id
                    print 'Processing ProteinID: ', g_id, '{:7,.0f}'.format(i)
                
                print_flag = True
                o = open(out+'/'+g_id +'.xml','w')
                o.write(string1)
                o.write('<interpro_matches>'+'\n')
            elif line.startswith('</protein>'):
                print_flag = False
                o.write(line+'\n')
                o.write('</interpro_matches>'+'\n')
                o.write(string2)
            if print_flag == True:
                o.write(line + '\n')

if __name__ == "__main__":
    
    options(sys.argv[1:])
    
    string1, string2 = IPRstring()
    
    readIn(string1, string2)
    
    ### close the logfile
    o.close()