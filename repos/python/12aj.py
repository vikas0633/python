### Usage: nice -n 19 python /Users/vikas0633/Desktop/plant/2012_week21/shortran_0.4/scripts/12aj.py /Users/vikas0633/Desktop/plant/2012_week21/shortran_0.4/demo "GHD-5_sample","GHD-6_sample","GHD-9_sample","GHD-10_sample" ouput_dir


### script for splitting the libraries

### importing modules							###
import time,os,datetime, sys,glob							###
work_dir=os.getcwd()							###
 
def split_one_folder(path,out):
	### make a data folder
	if not os.path.exists(out+'/module_1/'):
		os.mkdir(out+'/module_1/')
	if not os.path.exists(out+'/module_1/size_fractionated_data/'):
		os.mkdir(out+'/module_1/size_fractionated_data/')
	if not os.path.exists(out+'/module_1/size_fractionated_data/'+path.split('/')[-1]):
		os.mkdir(out+'/module_1/size_fractionated_data/'+path.split('/')[-1])
	
	
	### define a function for splitting files
	
	
	def split_files(file_name):
		### define an empty dictionary
		
		size_hash = {} ### for storing non-redundant size
		file = {} ### for storing file name variable
		
		### read fastq file and check for the size
		
		line_number = 0 ### tracker for the line number
		
		def write_file(temp_file,line1,line2,line3,line4):
			temp_file.write(line1+'\n')
			temp_file.write(line2+'\n')
			temp_file.write(line3+'\n')
			temp_file.write(line4+'\n')
		
		
		for line in open(file_name,'r'):
			line_number += 1
			line = line.strip()
			if (line_number %4 == 1):
				line1 = line
			elif (line_number %4 == 2):
				line2 = line
				size = len(line2)
			elif (line_number %4 == 3):
				line3 = line
			elif (line_number %4 == 0):
				line4 = line
			
				if size in size_hash:
					write_file(file[size],line1,line2,line3,line4)
				else:
					file[size]=open(out+'/module_1/size_fractionated_data/'+path.split('/')[-1]+'/'+str(size)+'.fastq','w')
					write_file(file[size],line1,line2,line3,line4)
					size_hash[size]=''
		return file
	### split the file
	for file_name in glob.glob( os.path.join(path, '*')):
		file = split_files(file_name)
	### close all the files
	for files in file:
		file[files].close()
		
		
### files to split
fq_path = sys.argv[1]
lib = sys.argv[2].split(',')
out = sys.argv[3]
for folder in lib:
	path=sys.argv[1]+'/'+folder
	split_one_folder(path,out)
		