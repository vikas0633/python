#12u.py - make a script which can make a fasta file with abundance in the header as of format Sequence_xAbundance

## for date in output file str(datetime.date.today())
import sys
import os
import glob
import re
import time
import datetime


### function for printing fasta file
def make_fasta_abundance(infile,column):
	line_count=0
	o=open(str(datetime.date.today())+'_'+(infile.split('/'))[-1],'w')
	for line in open(infile,'r'):
		line_count += 1
		if(len(line)>0):
			if line[0] in 'ATGCN':
				line=line.strip()
				token=line.split('\t')
				o.write(">"+str(line_count)+'_x'+str(token[int(column)-1])+'\n')
				o.write(token[0]+'\n')
	o.close()
if __name__ == "__main__":
	infile=raw_input("Enter the file name for making fasta file or press enter for default\n")
	print ("you have entered: ",infile)
	if(infile==''):
		infile='/Users/vikas0633/Desktop/plant/2011_week47/10_19_11_EINF_profile_19-24.sorted.cut-20.normalized'
		print ("Using input as: ",infile)
	column=raw_input("Enter the abundance column or press enter for default\n")
	print ("you have entered: ",column)
	if(column==''):
		column=2
		print ("Using input as: ",column)
	make_fasta_abundance(infile,column)