'''
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

### main script controling plant pipeline

### importing required modules
import os,sys,glob,re,commands,linecache,time,datetime

### ask for interactive configuration
print "Would you like to configure the run interactively ?"
print "press (Enter/Return) button for configering now"
print "press (n) if you have already changed the configuration file"
configure_run=raw_input("==>")

if (configure_run=="N" or configure_run=="n"):
	### remove .pyc file
	os.system('rm -rf `find . -name "*.pyc"`')
	### current directory as working directory
	work_dir=os.getcwd()
	lib_path = os.path.abspath(os.getcwd())
	script_dir=os.getcwd()+'/'+'scripts'
	sys.path.append(lib_path)
	import configuration_12
	date=configuration_12.date
	os.system('>time_stats.txt')
	#os.system('sh time_usage.sh >>time_stats.txt')

else:
	print "Please provide a folder name: "
	print "If default (today's date) press (enter/return)"
	print "If you would like to use previous analysis folder"
	print "provide the folder name without extention ['_output']"
	take=raw_input("==>")
	if (take==''):
		date=str(datetime.date.today()) 
	else:
		date=take
	print "These are the following modules can be run:"
	print "1.1 : Filtering and Sequencing Error correction"
	print "1.2 : making profiles"
	print "1.3 : Module for adding data and making mysql tables"
	print "2   : Genome mapping"
	print "3.1 : Cluster prediction"
	print "3.2 : Find genomic region for all the clusters"
	print "4.1 : Finding the miRNAs"
	print "4.2 : Map the miRNA to mirBase"
	print "5   : Ta-siRNA prediction"
	print "6   : Adding annotation to mysql database"
	print "7 : Making plot for expression values for libraries"
	print "8.1 : Find pattern in first three inucleotides, morkov probabilities"
	print "8.2 : Map the sequences with the pattern"
	print "9   : Make mySQL queries"
	
	print "Help: Provides a dumentation for modules"
	
	print "E/e : exit the program"
	print 
	print
	### remove .pyc file
	os.system('rm -rf `find . -name "*.pyc"`')
	### current directory as working directory
	work_dir=os.getcwd()
	lib_path = os.path.abspath(os.getcwd())
	script_dir=os.getcwd()+'/'+'scripts'
	sys.path.append(lib_path)
	import configuration_12i
	os.system('>time_stats.txt')
	#os.system('sh time_usage.sh >>time_stats.txt')
	
### funtion for module -1.1: Filtering and Sequencing Error correction
def Raw_data_filtering(configuration,output_dir):
	
	module_no=1
	### mapping information
	o=open(output_dir+'/filtering_summary_module_1','w')
	o.write("Input_file\t Number_of_reads\t Reasds_mapped_to_repeat\t Reads_not_mapped_to_repeat\n")
	### build index
	os.system('bowtie-build -q '+configuration.repeat_database+' '+ configuration.repeat_database +' 2> '+output_dir+'/module_1.logfile')			
	fq_file_dir = configuration.fq_file_dir
	### adapter filtering
	if (configuration.adapter_filtering == True):
		os.system('echo adapter filtering/trimming --- starts'+' >>time_stats.txt')
		t = time.time()
		adapter={};i=0
		if os.path.isfile(output_dir+'/module_1/adapter_filtering_summary.txt'):
				os.system('rm '+output_dir+'/module_1/adapter_filtering_summary.txt')
		if not os.path.exists(output_dir+'/module_1/'):
				os.mkdir(output_dir+'/module_1/')
		o_adapter=open(output_dir+'/module_1/adapter_filtering_summary.txt','a')
		o_adapter.write('File_name\t'+'Total_redundant_reads\t'+'Total_non-redundant_reads\t'+'Reads_contained_adapter_in_%\n')		
		for folder in configuration.lib:
			adapter[folder]=(configuration.adapter).split(',')[i]
			if len(configuration.adapter)==len(configuration.lib):
				i += 1
		for folder in configuration.lib:
			path=fq_file_dir+'/'+folder
			if not os.path.exists(output_dir+'/module_1/adapter_filtered_data/'):
				os.mkdir(output_dir+'/module_1/adapter_filtered_data/')
			if not os.path.exists(output_dir+'/module_1/adapter_filtered_data/'+path.split('/')[-1]):
				os.mkdir(output_dir+'/module_1/adapter_filtered_data/'+path.split('/')[-1])
			for infile in glob.glob( os.path.join(path, '*.f*q') ):
				if adapter[folder]!='':
					o_adapter=open(output_dir+'/module_1/adapter_filtering_summary.txt','a')
					line_count1=int((commands.getoutput('wc -l '+infile)).split()[0])
					line_count_unique=int((commands.getoutput('cat '+infile+' |sort|uniq -c| wc -l')).split()[0])
					#line_count=int(commands.getoutput('echo `cat '+miRNA_dataset+'|wc -l`'))
						
					os.system('chmod 777 '+ script_dir +'/*')
					### check if mac
						
					if sys.platform=='darwin':	
						if configuration.keep_with_adapter_only==True:
							os.system('nice -n 19 '+ script_dir +'/fastx_clipper_osx -Q 33 -c -a '+adapter[folder]+' -i '+infile + ' -o '+output_dir+'/module_1/adapter_filtered_data/'+path.split('/')[-1]+'/'+infile.split('/')[-1])	
						else:
							os.system('nice -n 19 '+ script_dir +'/fastx_clipper_osx -Q 33 -a '+adapter[folder]+' -i '+infile + ' -o '+output_dir+'/module_1/adapter_filtered_data/'+path.split('/')[-1]+'/'+infile.split('/')[-1])
					else:
					
						if configuration.keep_with_adapter_only==True:
							os.system('nice -n 19 '+ script_dir +'/fastx_clipper -Q 33 -c -a '+adapter[folder]+' -i '+infile + ' -o '+output_dir+'/module_1/adapter_filtered_data/'+path.split('/')[-1]+'/'+infile.split('/')[-1])	
						else:
							os.system('nice -n 19 '+ script_dir +'/fastx_clipper -Q 33 -a '+adapter[folder]+' -i '+infile + ' -o '+output_dir+'/module_1/adapter_filtered_data/'+path.split('/')[-1]+'/'+infile.split('/')[-1])
					infile=output_dir+'/module_1/adapter_filtered_data/'+path.split('/')[-1]+'/'+infile.split('/')[-1]	
					line_count2=int((commands.getoutput('wc -l '+infile)).split()[0])
					o_adapter.write(infile.split('/')[-1]+'\t'+str(line_count1/4)+'\t'+str(line_count_unique/4)+'\t'+str(round(float(line_count2)*100/float(line_count1),2))+'\n')
				else:
					os.system('cp '+infile+' '+output_dir+'/module_1/adapter_filtered_data/'+path.split('/')[-1]+'/'+infile.split('/')[-1])
					infile=output_dir+'/module_1/adapter_filtered_data/'+path.split('/')[-1]+'/'+infile.split('/')[-1]
				
				#########################
				### Trimming ############
				#########################
				
				trim = int(configuration.trim)
				if trim != 0:
					if sys.platform=='darwin':
						os.system('nice -n 19 '+ script_dir +'/fastx_trimmer_osx -Q 33 -t '+str(trim)+' -i '+infile + ' -o '+infile+'_trimmed_'+str(trim)+'.fastq')
					else:	
						os.system('nice -n 19 '+ script_dir +'/fastx_trimmer -Q 33 -t '+str(trim)+' -i '+infile + ' -o '+infile+'_trimmed_'+str(trim)+'.fastq')
					
					
					os.system('rm '+infile)	
			o_adapter.close()
		fq_file_dir = output_dir+'/module_1/adapter_filtered_data/'
		
		
		### plot the abundances
		os.system('nice -n 19 python '+script_dir+'/12am.py '+output_dir+'/module_1/adapter_filtering_summary.txt')
		os.system('mv *.png '+output_dir)
		os.system('echo adapter filtering/trimming --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')	
	#print resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
	### fastq size fractionation
	if (configuration.size_fractionation == True):
		os.system('echo size splitting --- starts'+' >>time_stats.txt')
		t = time.time()
		lib = ''
		for l in configuration.lib:
			lib += l +','
		lib = lib[:-1]
		os.system('nice -n 19 python '+script_dir+'/12aj.py '+fq_file_dir+' '+lib+' '+output_dir)
		fq_file_dir = output_dir+'/module_1/size_fractionated_data'
		os.system('echo size splitting --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')
	
	t = time.time()
	os.system('echo filtering/error correction --- starts'+' >>time_stats.txt')	
	for folder in configuration.lib:
		path=fq_file_dir+'/'+folder
		for infile in glob.glob( os.path.join(path, '*.*f*q') ):
			### move the files from .fq to .fastq version for Echo
			import re
			if re.search('.fq',infile):
				outfile = infile.replace('.fq','.fastq')
				os.system('mv '+infile+' '+ outfile)
				infile = outfile
		#for infile in glob.glob( os.path.join(path, '*corrected') ):
			#o.write("Input file: ")
			### if file name consist correct size
			size=int(commands.getoutput('echo `awk ''NR==2'' '+infile+' |wc -c`'))
			size=size-1
			#size=len((linecache.getline(infile, 2)).strip())
			if ( configuration.min_seq_size <= int(size) <= configuration.max_seq_size):
			
				if(configuration.ErrorCorrection==True):
					os.system('nice -n 19 python '+ configuration.Path_2_ECHO +'/ErrorCorrection.py -o ' + infile+'_corrected ' +infile+' 1> '+output_dir+'/module_1.logfile'+' 2> '+output_dir+'/module_1.logfile') 
					os.system('rm '+path+'/*.qual')
					os.system('rm '+path+'/*.seq')
					infile = infile+'_corrected'

			
				o.write(infile.split('/')[-1])
				o.write('\t')
				### map the reads
				os.system('nice -n 19 bowtie --threads '+str(configuration.no_of_cores)+' --seedmms '+str(configuration.mismatches)+' --maqerr '+str(configuration.max_maq_err)+' --seedlen '+str(configuration.seedlen)+' -k 1 --offbase 0 --sam '+configuration.repeat_database+' '+ infile+' '+' 2>> '+output_dir+'/module_1.logfile'+" | nice -n 19 perl -lane 'print $_ if ($F[3] > 0)' >"+ output_dir+'/'+(infile.split('/'))[-1]+'.sam ')
				### remove reads which were mapped on repeats
				os.system('nice -n 19 python '+script_dir+'/12m.py'+' '+output_dir+'/'+(infile.split('/'))[-1]+'.sam'+' '+infile +' '+output_dir+'/'+folder+'_'+(infile.split('/'))[-1]+"_repeat_filtered > "+output_dir+'/'+folder+'_'+(infile.split('/'))[-1]+"_repeat")
				
				#o.write("number of sequences in the input fastq file: ")
				number=sum(1 for line in open(infile))
				o.write(str(int(number)/4))
				o.write('\t')
				#o.write("number of sequences mapped to repeat: ")
				number=sum(1 for line in open(output_dir+'/'+folder+'_'+(infile.split('/'))[-1]+"_repeat"))
				o.write(str(int(number)/4))
				o.write('\t')
				#o.write("number of sequences not mapped to repeat: ")
				number=sum(1 for line in open(output_dir+'/'+folder+'_'+(infile.split('/'))[-1]+"_repeat_filtered"))
				o.write(str(int(number)/4))
				o.write('\n')
	o.close()
	os.system('rm '+output_dir+'/'+'*sam')
	os.system('echo filtering/error correction --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')
	if(configuration.use_only_mapped_reads==True):
		os.system('rm '+output_dir+'/'+'*repeat_filtered')
	else:
		os.system('rm '+output_dir+'/'+'*repeat')
	
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
		
	os.system('mv '+output_dir+'/*module_1* '+output_dir+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
	os.system('mv '+output_dir+'/*repeat* '+output_dir+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')

### funtion for module -1.2: making profiles
def Make_profile_files(configuration,output_dir):
	import os
	t = time.time()
	os.system('echo making profile --- starts >>time_stats.txt')
	module_no=1
	#print "libraries: ",configuration_12.lib
	### make a loop for all size classes and make output file
	os.system('>'+output_dir+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff))
	for i in range(configuration.min_seq_size, configuration.max_seq_size+1,1):
		os.system('echo "#Sequence" >'+output_dir+'/'+date+'_profile_'+str(i))
		for folder in configuration.lib:
			path=output_dir+'/module_'+str(module_no)
			if(configuration.use_only_mapped_reads==True):
				for infile in glob.glob( os.path.join(path, '*'+folder+'*repeat')):
					size=int(commands.getoutput('echo `awk ''NR==2'' '+infile+' |wc -c`'))
					size=size-1
					if (size==i):
			
						### making profile using python
						os.system('nice -n 19 python '+script_dir+'/12l.py '+infile+' '+output_dir+'/'+ date+'_profile_'+str(i))
				
			
			else:
			
				for infile in glob.glob( os.path.join(path, '*'+folder+'*repeat_filtered')):
					size=int(commands.getoutput('echo `awk ''NR==2'' '+infile+' |wc -c`'))
					size=size-1
					if (size==i):
						### making profile using python
						os.system('nice -n 19 python '+script_dir+'/12l.py '+infile+' '+output_dir+'/'+ date+'_profile_'+str(i))
				
		### sort the profiles
		#print "awk 'NR == 1; NR > 1 {print $0 | "+'"sort"}'+"' "+ str(datetime.date.today())+'_profile_'+str(i)+ ' > '+str(datetime.date.today())+'_profile_'+str(i)+'_sorted'
		os.system("awk 'NR == 1; NR > 1 {print $0 | "+'"sort"}'+"' "+output_dir+'/'+ date+'_profile_'+str(i)+ ' > '+output_dir+'/'+date+'_profile_'+str(i)+'_sorted')
		
		
		### take the header
		os.system('head -1 '+output_dir+'/'+date+'_profile_'+str(i)+'_sorted'+' > '+output_dir+'/'+date+'_profile_'+str(i)+'_sorted'+'.cut-'+str(configuration.cutoff))
		### count abundances in each column and print to file called check
		os.system("awk '{ s = 0; for (i = 1; i <= NF; i++) s = s+$i;"+ 'print $0"\t"s }'+"' "+output_dir+'/'+date+'_profile_'+str(i)+'_sorted'+ '> check')
		os.system('awk -f '+script_dir+'/sum_all.awk check >'+output_dir+'/'+ 'lib_abundance_'+str(i))
		os.system("awk '{ s = 0; for (i = 1; i <= NF; i++) s = s+$i; if (s > "+str(configuration.cutoff)+') print $0"\t"s }'+"' "+output_dir+'/'+ date+'_profile_'+str(i)+'_sorted'+'>>'+output_dir+'/'+ date+'_profile_'+str(i)+'_sorted'+'.cut-'+str(configuration.cutoff))
		os.system('rm check')
		
		if(i==configuration.min_seq_size):
			os.system('cat '+output_dir+'/'+ date+'_profile_'+str(i)+'_sorted'+'.cut-'+str(configuration.cutoff)+' >> '+output_dir+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff))
			os.system('nice -n 19 python '+script_dir+'/12ab.py '+'>'+output_dir+'/libraries_abundance')
		else:
			os.system('sed 1d '+output_dir+'/'+ date+'_profile_'+str(i)+'_sorted'+'.cut-'+str(configuration.cutoff)+' >> '+output_dir+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff))
		
		### make size distribution table that will be used for making some nice pie charts
		os.system('nice -n 19 python '+script_dir+'/12ac.py '+output_dir+'/'+ 'lib_abundance_'+str(i) +' '+str(i)+' >>'+output_dir+'/libraries_abundance')

		print "processed size "+str(i)
				
		os.remove(output_dir+'/'+ 'lib_abundance_'+str(i))
		## remove files
		os.system('rm '+output_dir+'/'+ date+'_profile_'+str(i))
		os.system('rm '+output_dir+'/'+date+'_profile_'+str(i)+'_sorted')

			
	
	### normalization using python
	infile = output_dir+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)
	os.system('nice -n 19 python '+script_dir+'/12n.py '+infile+ '>'+infile+'_normalized')
	### make files for each size
	for i in range(configuration.min_seq_size, configuration.max_seq_size+1,1):
		os.system('head -1 '+infile+'_normalized' +' > '+output_dir+'/'+ date+'_profile_'+str(i)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized')
		os.system(" awk ' {if(length($1)=="+str(i)+") print $0} ' "+ infile+'_normalized' +' >> '+output_dir+'/'+ date+'_profile_'+str(i)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized')
	infile = infile+'_normalized'
	outfile=infile+'.fasta'
	## check if infile is empty
	number=sum(1 for line in open(infile))
	if number == 0:
		sys.exit('profile table is empty')
		
	
	### making fasta files with the sequences
	os.system('nice -n 19 python '+script_dir+'/12c.py '+ infile +' '+ outfile)
	
	### apply additional score cut-off
	os.system('nice -n 19 python '+script_dir+'/12w.py '+ infile+' '+str(configuration.score_regulation_cutoff)+' > '+infile+'_score_cutoff_'+str(configuration.score_regulation_cutoff))
	
	### make a file with unique count file
	os.system('nice -n 19 python '+ script_dir+'/12ai2.py '+ infile+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+' '+output_dir+'/libraries_abundance')
	os.system('nice -n 19 python '+ script_dir+'/12ai.py '+ infile+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+' '+output_dir+'/libraries_abundance_unique')
	
	### plot the abundances
	os.system('nice -n 19 python '+script_dir+'/12ad.py '+output_dir+'/libraries_abundance'+' '+output_dir+'/libraries_abundance_unique')

	### make compare library files
	import os.path
	if os.path.isfile(output_dir+'/'+'libraries_reference'):
		print 'library references are provided'
	else:
		os.system('nice -n 19 python '+script_dir+'/12o.py '+str(len(configuration.lib))+' '+output_dir+'/'+'libraries_reference > '+output_dir+'/'+'compare_libraries')

	### move plots
	if not os.path.exists(output_dir+'/module_'+str(module_no)+'/plots'):
		os.mkdir(output_dir+'/module_'+str(module_no)+'/plots')
	os.system('mv '+str(date)+'_output/*png '+str(date)+'_output'+'/module_'+str(module_no)+'/plots' + ' 2>'+output_dir+'/STDOUT.txt')
		
	### move profiles
	if not os.path.exists(output_dir+'/module_'+str(module_no)+'/profiles'):
		os.mkdir(output_dir+'/module_'+str(module_no)+'/profiles')
	os.system('mv '+str(date)+'_output/*profile* '+str(date)+'_output'+'/module_'+str(module_no)+'/profiles' + ' 2>'+output_dir+'/STDOUT.txt')

	
	import os
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
			
	if not os.path.exists(output_dir+'/module_'+str(module_no)+'/filter'):
		os.mkdir(output_dir+'/module_'+str(module_no)+'/filter')
	os.system('mv '+output_dir+'/module_'+str(module_no)+'/*repeat* '+output_dir+'/module_'+str(module_no)+'/filter/' )
	os.system('mv '+output_dir+'/*libraries* '+output_dir+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
		
	### python script for reaplcing libraries header
	infile=output_dir+'/module_'+str(module_no)+'/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
	os.system('cp '+infile+' '+'check.txt')
	os.system('nice -n 19 python '+script_dir+'/12ae.py '+'check.txt '+infile)
	
	os.system('echo making profile --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')

### funtion for module -1.3: module for adding data and making mysql tables
def add_output_to_mysql_database(configuration,output_dir):
	import os
	t = time.time()
	os.system('echo adding to mysql --- starts >>time_stats.txt')
	module_no=1
	if(len(configuration.infile)==0):
		infile=output_dir+'/module_'+str(module_no)+'/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
	else:
		infile=configuration.infile
	### add the length column
	o=open(infile+'_len','w')
	first_line=False
	for line in open(infile,'r'):
		line=line.strip()
		token=line.split('\t')
		if first_line==False:	
			string=token[0]+'\t'
			string += 'read_size\t'
		else:
			string=token[0]+'\t'
			string += str(len(token[0]))+'\t'
		for i in range(1,len(token)):
			string += token[i]+'\t'
		o.write(string+'\n')
		first_line=True
	o.close()
	os.system('sed 1d '+infile+'_len'+' >check.txt')
	os.system('nice -n 19 python '+script_dir+'/12r.py '+infile+'_len' +' '+work_dir+' '+output_dir+'/temp_folder/index.txt >'+output_dir+'/temp_folder/mysql_batch')
	os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/mysql_batch')
	os.system('rm check.txt')
	os.system( 'rm '+infile+'_len')
	os.system('echo adding to mysql --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')

### funtion for module -2: Genome mapping
def genome_mapping(configuration,output_dir):
	import os
	t=time.time()
	os.system('echo genomic --- starts >>time_stats.txt')
	module_no=2
	if(len(configuration.infile)==0):
		infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
	else:
		infile=configuration.infile
	file=configuration.genome_file
	
	outfile=infile+'.fasta'
	os.system('nice -n 19 python '+script_dir+'/12c.py '+ infile +' '+ outfile)
	os.system("bowtie-build "+file+' '+file+' 1>'+output_dir+'/map.logfile')
	os.system('nice -n 19 bowtie --threads '+str(configuration.no_of_cores)+' --seedmms '+str(configuration.mismatches)+' --maqerr '+str(configuration.max_maq_err)+' --seedlen '+str(configuration.seedlen)+' -k 30 --best --time --offbase 0 --sam '+ file +' -f '+ outfile +' 2>'+output_dir+'/map.logfile'+"| nice -n 19 perl -lane '"+"print $_ if ($F[3] > 0)' > "+ infile+'.sam ')
	
	### infile.sam will have all the mapping positions on genome
	### we need to combine normalized expression file with genome mappings
	os.system('nice -n 19 python '+script_dir+'/12v.py '+infile+' '+outfile+' '+infile+'.sam > '+infile+'.igv')
	os.system('head -n 1 '+infile+'.igv > '+infile+'.sorted.igv')
	os.system('nice -n 19 sed 1d '+infile+'.igv '+'|sort -b -k 1,1  -k 2,2n >> '+infile+'.sorted.igv')
	os.system('rm '+infile+'.igv')
	
	### check if IGV file is empty	
	number=sum(1 for line in open(infile+'.sorted.igv'))
	if number == 0:
		sys.exit('Mapped IGV is empty')

	os.system('echo Genomic mapping --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')
	### move files to module-2 folder
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	os.system('mv '+outfile+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
	os.system('mv '+infile+'.sorted.igv'+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
	os.system('mv '+output_dir+'/map.logfile'+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
	os.system('mv '+infile+'.sam'+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')


### funtion for module -3.1: cluster prediction
def cluster_prediction(configuration,output_dir):
	module_no=3
	import os
	t=time.time()
	os.system('echo cluster prediction --- starts >>time_stats.txt')
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))

	if(len(configuration.igv_infile)==0):
		infile=output_dir+'/module_2'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+'.sorted.igv'
	else:
		infile=configuration.igv_infile
	output=output_dir+'/module_2'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+'.sorted.igv'+'_cluster.bed'
	
	### cluster prediction
	if configuration.NiBLS == True:
		os.system('nice -n 19 perl '+script_dir+'/nibls.pl -f ' +infile +'>'+output+'.gff '+' 2>>'+output_dir+'/ERROR_OUT.txt')
		os.system('cut -f 1,4-6 '+output+'.gff >'+output)
	else:
		os.system('nice -n 19 python '+script_dir+'/12g.py ' +infile +' '+output+' '+str(configuration.divide_by_size)+' '+str(configuration.abundance_cutoff)+' '+str(configuration.min_seq_size)+' '+str(configuration.max_seq_size)+' '+str(configuration.max_gap_size)+' '+str(configuration.min_unique_reads)+' > cluster_summary.txt'+' 2>>'+output_dir+'/ERROR_OUT.txt')
	
	### check if cluster Bed file is empty	
	number=sum(1 for line in open(output))
	if number == 0:
		sys.exit('cluster Bed file is empty')
		
	os.system('mv '+output_dir+'/module_2/*.bed*'+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>>'+output_dir+'/STDOUT.txt')
		
	
	if(configuration.Add_cluster_to_mySQL==True):
		infile = output_dir+'/module_3'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+'.sorted.igv'+'_cluster.bed'
		#--#os.system('cp '+infile+' check.txt')
		#--#os.system('nice -n 19 python '+script_dir+'/12r_cluster.py '+infile +' '+work_dir+' '+output_dir+'/temp_folder/index_cluster.txt >'+output_dir+'/temp_folder/mysql_batch')
		#--#os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/mysql_batch')
		### make cluster fasta file
		### Make blast database for genome
		output=output_dir+'/module_3'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+'.sorted.igv'+'_cluster.bed'
		os.system('cp '+output+' '+'clusters')
		output='clusters'
		os.system('>'+output+'.fa')
		os.system('nice -n 19 formatdb -i '+configuration.genome_file+' -o T -p F')
		
		
		### open a file to store the coordinates for format db
		fdb=open('Cluster_loci.sh','w')
		first_line=True
		for line in open(output,'r'):
			line=line.strip()
			if first_line==False:			
				token=line.split("\t")
				### writing cluster sequences
				fdb.write('fastacmd -d '+configuration.genome_file+' -p F -s '+token[0]+" -L "+token[1]+','+token[2]+'\n')
			first_line=False
		fdb.close()
		os.system('sh Cluster_loci.sh'+' >>'+output+'.fa'+' 2>'+output_dir+'/STDOUT.txt')
		os.system('rm Cluster_loci.sh')
		if(len(configuration.infile)==0):
			infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
		else:
			infile=configuration.infile
		outfile=infile+'.fasta'
		os.system('nice -n 19 python '+script_dir+'/12c.py '+ infile +' '+ outfile+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt')
		file=output+'.fa'
		os.system('>'+output_dir+'/temp_folder/mysql_batch')
		### index the reference
		os.system("bowtie-build -q "+file+' '+file+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt')
		### map the reads
		os.system('nice -n 19 bowtie -q --threads '+str(configuration.no_of_cores)+' --seedmms '+str(configuration.mismatches)+' --maqerr '+str(configuration.max_maq_err)+' --seedlen '+str(configuration.seedlen)+' -k 100 --best --time --offbase 0 --sam '+ file +' -f '+ outfile+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt' +"| nice -n 19 perl -lane '"+"print $_ if ($F[3] > 0)' "+'|sort -nr -k 5,5 '+" > "+ infile+'.sam')
		### make file with sequences and annotations
		os.system('nice -n 19 python '+script_dir+'/12s.py '+infile+' '+ infile+'.sam '+ file + ' >'+output_dir+'/'+'out.txt'+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt')
		os.system('nice -n 19 python '+script_dir+'/12t.py '+infile+' '+output_dir+'/'+'out.txt '+ file +' >> '+output_dir+'/temp_folder/mysql_batch'+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt')
		os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/mysql_batch'+'\n')
		os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/index.txt'+'\n')
	
		
	#os.system('rm clusters*')
	#os.system('rm clusters.fa')
	### move files to module-3 folder
	os.system('mv '+output+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
	os.system('mv cluster_summary.txt'+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
	os.system('mv '+output_dir+'/ERROR_OUT.txt'+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
	os.system('mv '+infile+'.sam'+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
	os.system('echo cluster prediction --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')

### funtion for module -3.X: calculating unique sequences in the clusters
def calculate_unique_sequences_in_cluster(configuration,output_dir):
	module_no=3
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	if(len(configuration.cluster_infile)==0):
		infile=output_dir+'/module_3'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+'.sorted.igv_cluster.bed'
	else:
		infile=configuration.cluster_infile
	o=open(output_dir+'/module_3/unique_sequence_counts.txt','w')
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	o.write("number of clusters: "+ str(int(commands.getoutput('echo `cat '+infile+'|wc -l`'))-1)+'\n')
	o.close()
	os.system('nice -n 19 python '+script_dir+'/12k.py '+ infile +'>> '+output_dir+'/module_3/unique_sequence_counts.txt' ) 

### funtion for module -3.X: make a file with regulated sequences and find these sequences in clusters
def regulated_sequences_in_cluster(configuration,output_dir):
	module_no=3
	### check if regulation file is required
	if(len(configuration.regulation_file)==0):
		infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'
		os.system('nice -n 19 python '+script_dir+'/12w.py '+ infile+' '+str(configuration.score_regulation_cutoff)+' > '+infile+'_score_cutoff_'+str(configuration.score_regulation_cutoff))
	else:
		infile=configuration.regulation_file
	### check if cluster file is provided
	if(len(configuration.cluster_infile)==0):
		cluster_infile=output_dir+'/module_3'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+'.sorted.igv_cluster.bed'
	else:
		cluster_infile=configuration.cluster_infile
	os.system('nice -n 19 python '+script_dir+'/12i.py '+infile+' '+cluster_infile +'>> '+output_dir+'/module_3/regulated_sequence_counts.txt')

### funtion for module -3.2: find genomic region for all the clusters
def cluster_genomic_region_prediction(configuration,output_dir):
	import os
	t=time.time()
	os.system('echo cluster_genomic_region --- starts >>time_stats.txt')
	module_no=3
	if(len(configuration.cluster_infile)==0):
		cluster_infile=output_dir+'/module_3'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+'.sorted.igv_cluster.bed'
	else:
		cluster_infile=configuration.cluster_infile
	os.system('nice -n 19 python '+script_dir+'/12h.py '+configuration.gtf_file+' '+cluster_infile +' > '+cluster_infile+'.region')
	exonic=0;intronic=0;overlap=0;intergenic=0;
	for line in open(cluster_infile+'.region','r'):
		line = line.strip()
		token=line.split('\t')
		if (token[5]=='exonic'):
			exonic += 1
		if (token[5]=='intronic'):
			intronic += 1
		if (token[5]== "exon-intron overlap"):
			overlap +=1
		if (token[5]== "intergenic"):
			intergenic +=1
			
	sum=exonic+intronic+intergenic+overlap
		
	### making pichart for cluster genomic distribution
	import matplotlib
	matplotlib.use('Agg')
	import matplotlib.pyplot as plt
	plt.pie([exonic/float(sum),intronic/float(sum),intergenic/float(sum),overlap/float(sum)],colors=((0.2023,1,0.8),(0.2023,0.8,1),(1,0.4745,0.3),(0.8509,0.4,1)),shadow=True,labels=['exonic','intronic','intergenic','overlapping'],pctdistance=0.6, labeldistance=1.1,autopct='%1.f%%',explode=(0.0,0.1,0.1,0.05))
	plt.title('Genomic Distribution of Clusters')
	plt.savefig(cluster_infile+'.region.png')

	o=open('cluster_genomic_region.txt','w')
	o.write("number of exonic cluster: "+str(exonic)+'\n')
	o.write("number of intronic cluster: "+str(intronic)+'\n')
	o.write("number of intergenic cluster: "+str(intergenic)+'\n')
	o.write("number of exon-intron overlapping cluster: "+str(overlap)+'\n')
	o.close()
	
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	os.system('mv '+'cluster_genomic_region.txt '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')

	os.system('echo cluster_genomic_region --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')
	
### funtion for module -4.1:	finding the miRNAs
def miRNA_prediction(configuration,output_dir):
	import os
	t=time.time()
	os.system('echo miRNA predictions --- starts >>time_stats.txt')
	module_no=4
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	dir=configuration.miRNA_dir
	mi_genome_file=configuration.mi_genome_file
	if(len(configuration.mi_genome_file)==0):
		print "please provide correct reference file for miRNA predictions"
	else:
		###index the genome
		os.system('cp '+mi_genome_file+' '+mi_genome_file.split('/')[-1])
		#print 'bowtie-build -f '+mi_genome_file+' '+ mi_genome_file
		os.system('nice -n 19 bowtie-build -q -f '+mi_genome_file.split('/')[-1]+' '+ mi_genome_file.split('/')[-1] +' 2>'+output_dir+'/miRNA_map.txt')		
		if configuration.mirDeep_2==True:
			### making miRNA predictions using miRDeep 2.0
			if(len(configuration.mi_sequence_file)==0):
				o=open(output_dir+'/module_'+str(module_no)+'/miRNA','w')
				infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
				c=0
				for line in open(infile,'r'):
					line=line.strip()
					if(len(line)>0):
						c += 1
						token=line.split('\t')
						if(line[0]=='#'):
							i = token.index('sum')
						else:
							###>ShootPi1000000_x3
							o.write('>'+'seq'+'_'+str(c)+'_x'+str((int(token[i])/len(configuration.lib))+1)+'\n')
							o.write(token[0]+'\n')
				o.close()
				infile=output_dir+'/module_'+str(module_no)+'/miRNA'
				
				### only for c. elegens
				#infile = "/Users/vikas0633/Desktop/plant/2012_week22/mirdeep2/tutorial_dir/reads_collapsed.fa"
				
				
			else:
				infile=mi_sequence_file

			
			
			### use mapper.pl to make mapped .arf file
			os.system('cp '+infile+' miRNA')
			os.system('mapper.pl '+'miRNA'+' -c -p '+ mi_genome_file.split('/')[-1] +' -t '+'miRNA'+'.arf ' +' >> '+output_dir+'/miRNA_map.txt')
			### use mirDeep 2.0 to do the prediction
			os.system('miRDeep2.pl '+ 'miRNA' +' '+mi_genome_file.split('/')[-1]+' '+'miRNA'+'.arf '+configuration.miR_known_mature+' '+configuration.miR_related_mature+' '+configuration.miR_known_precursor+' 2> report.log')
			
			## make output to parse file in mysql
			os.system('nice -n 19 python '+script_dir+'/12al.py -i result_*.csv -o '+output_dir+'/module_'+str(module_no)+'/miRNA_predictions 2>>'+output_dir+'/STDOUT.txt')
		
			### making output that can be added to the MySQL table 	
			os.system('>'+output_dir+'/temp_folder/mysql_batch')
			infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
			os.system('nice -n 19 python '+script_dir+'/12af.py '+infile+' '+output_dir+'/module_'+str(module_no)+'/miRNA_predictions'+' '+output_dir+'/'+'check.txt '+' 2>>'+output_dir+'/STDOUT.txt')
			os.system('nice -n 19 python '+script_dir+'/12t.py '+infile+' '+output_dir+'/'+'check.txt '+ 'miRNA_prediction_mirdeep2' +' >>'+output_dir+'/temp_folder/mysql_batch'+' 2>>'+output_dir+'/STDOUT.txt')
			os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/mysql_batch'+' 2>>'+output_dir+'/STDOUT.txt')
			os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/index.txt'+' 2>'+output_dir+'/STDOUT.txt')
	
			### move .txt file to module-4 folder
			os.system('mv '+output_dir+'/*.txt '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>>'+output_dir+'/STDOUT.txt')	
			os.system('rm '+output_dir+'/module_'+str(module_no)+'/miRNA'+' 2>>'+output_dir+'/STDOUT.txt')	
			os.system('mv mirdeep_runs '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>>'+output_dir+'/STDOUT.txt')
			os.system('mv expression_analyses '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>>'+output_dir+'/STDOUT.txt')
			os.system('mv dir* pdfs* '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>>'+output_dir+'/STDOUT.txt')	
			os.system('mv *.log *.html *.arf *.csv '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>>'+output_dir+'/STDOUT.txt')
			os.system('rm *.ebwt')
		else:
		
			if(len(configuration.mi_sequence_file)==0):
				o=open(output_dir+'/module_'+str(module_no)+'/miRNA','w')
				infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
				for line in open(infile,'r'):
					line=line.strip()
					if(len(line)>0):
						token=line.split('\t')
						if(line[0]=='#'):
							i = token.index('sum')
						else:
							###>ShootPi1000000_x3
							o.write('>'+token[0]+'_x'+str((int(token[i])/len(configuration.lib))+1)+'\n')
							o.write(token[0]+'\n')
				o.close()
				infile=output_dir+'/module_'+str(module_no)+'/miRNA'
			else:
				infile=mi_sequence_file

		
		
			### make chromosome_length.txt file
			os.system('nice -n 19 python '+script_dir+'/12x.py ' + mi_genome_file +' >'+output_dir+'/chromosome_length.txt')
			
			#print int(commands.getoutput('echo `cat '+infile+'|wc -l`'))
			### map the reads
			#print 'bowtie -v 0 '+mi_genome_file+' -f '+ infile+' > '+infile+'.aln'
			os.system('nice -n 19 bowtie --threads '+str(configuration.no_of_cores)+' -v 0 '+mi_genome_file+' -f '+ infile+' > '+infile+'.aln' +' 2>>'+output_dir+'/miRNA_map.txt')
			
			###convert aln format into blast format
			#print 'perl '+dir+'/convert_bowtie_to_blast.pl '+infile+'.aln' +' '+ infile +" "+mi_genome_file+' > '+infile+'.bst'
			os.system('nice -n 19 perl '+dir+'/convert_bowtie_to_blast.pl '+infile+'.aln' +' '+ infile +" "+mi_genome_file+' > '+infile+'.bst')
	
			###FILTERING UBIQUITOUS ALIGNMENT
			#print 'perl '+dir+'/filter_alignments.pl '+infile+'.bst'+' -c 15 > '+ infile+'_filtered15.bst'
			os.system('nice -n 19 perl '+dir+'/filter_alignments.pl '+infile+'.bst'+' -c 15 > '+ infile+'_filtered15.bst')
			###For other species, different cutoffs, based on the known family sizes or other empirical 
			###considerations such as genome sizes, might be selected
			
			###FILTERING ALIGNMENTS BY ANNOTATION
			###The reads which mapped to exons or other non-coding RNAs, including rRNA, snRNA, snoRNA, 
			###and tRNA, are filtered by:
			if(len(configuration.gff_file)==0):
				os.system('>'+output_dir+'/gff_test.txt')
				configuration.gff_file=output_dir+'/gff_test.txt'
				print 'did you forget to provide annotation file ??  '
			#print 'perl '+dir+'/overlap.pl '+infile+'_filtered15.bst' +' '+configuration.gff_file+' -b >'+ infile+'_id_overlap_ncRNA_CDS'
			os.system('nice -n 19 perl '+dir+'/overlap.pl '+infile+'_filtered15.bst' +' '+configuration.gff_file+' -b >'+ infile+'_id_overlap_ncRNA_CDS') 
			
			###Only alignments where the read ids are not included in the id_overlap_ncRNA_CDS file are 
			###retained. g designates that lines where the query read ids are included in the 
			###id_overlap_ncRNA_CDS file are discarded
			#print 'perl '+dir+'/alignedselected.pl '+infile+'_filtered15.bst'+' -g '+infile+'_id_overlap_ncRNA_CDS > '+infile+'_filter15_ncRNA_CDS.bst'
			os.system('nice -n 19 perl '+dir+'/alignedselected.pl '+infile+'_filtered15.bst'+' -g '+infile+'_id_overlap_ncRNA_CDS > '+infile+'_filter15_ncRNA_CDS.bst')
			###Reads can also be filtered such that only reads that have one or more remaining alignments are 
			###kept. -b designates that the output should be FASTA entries and not alignments:
			#print 'perl '+dir+'/filter_alignments.pl '+infile+'_filter15_ncRNA_CDS.bst'+' -b '+infile+' > '+infile+'_filtered.fa'
			os.system('nice -n 19 perl '+dir+'/filter_alignments.pl '+infile+'_filter15_ncRNA_CDS.bst'+' -b '+infile+' > '+infile+'_filtered.fa')
			
			### PREDICTING MICRORNAS
			###Using the remaining alignments as guidelines, the potential precursor sequences are excised 
			###from the genome. This step is time-consuming, especially when the reference genome is large.
			#print 'perl '+dir+'/excise_candidate.pl '+mi_genome_file+' '+  infile+'_filter15_ncRNA_CDS.bst'+' 250 >'+infile+'_precursors.fa'
			os.system('nice -n 19 perl '+dir+'/excise_candidate.pl '+mi_genome_file+' '+  infile+'_filter15_ncRNA_CDS.bst'+' 250 >'+infile+'_precursors.fa')
			
			###The secondary structures of the potential precursors are predicted using RNAfold. -noPS means that no graphical output is produced. 
			#print 'cat '+infile+'_precursors.fa'+' | RNAfold -noPS > '+infile+'_structures'
			os.system('nice -n 19 cat '+infile+'_precursors.fa'+' | RNAfold -noPS > '+infile+'_structures')
			
			###The signatures are generated by aligning the remaining reads to the potential precursors.
			#print 'bowtie-build -f '+infile+'_precursors.fa '+ infile+'_precursors.fa'
			os.system('nice -n 19 bowtie-build -q -f '+infile+'_precursors.fa '+ infile+'_precursors.fa'+' 2>>'+output_dir+'/miRNA_map.txt') 
			#print 'bowtie -a -v 0 '+infile+'_precursors.fa '+' -f '+infile+'_filtered.fa'+' > '+infile+'_precursors.aln' 
			os.system('nice -n 19 bowtie -a -v 0 --threads '+str(configuration.no_of_cores)+' '+infile+'_precursors.fa '+' -f '+infile+'_filtered.fa'+' > '+infile+'_precursors.aln' +' 2>>'+output_dir+'/miRNA_map.txt')
			#print 'perl '+dir+'/convert_bowtie_to_blast.pl '+infile+'_precursors.aln'+' '+ infile+'_filtered.fa '+ infile+'_precursors.fa'+' > '+infile+'_precursors.bst'
			os.system('nice -n 19 perl '+dir+'/convert_bowtie_to_blast.pl '+infile+'_precursors.aln'+' '+ infile+'_filtered.fa '+ infile+'_precursors.fa'+' > '+infile+'_precursors.bst') 
			
			###Note that it is necessary for the prediction to do this sorting. 
			#print 'sort +3 -25 '+infile+'_precursors.bst'+' > '+infile+'_signatures'
			os.system('nice -n 19 sort -k3,25 '+infile+'_precursors.bst'+' > '+infile+'_signatures') 
			
			###Predictions are made
			#print 'perl '+dir+'/miRDP.pl '+infile+'_signatures '+ infile+'_structures'+' -y > '+infile+'_predictions'
			os.system('nice -n 19 perl '+dir+'/miRDP.pl '+infile+'_signatures '+ infile+'_structures'+' -y > '+infile+'_predictions') 
			
			###REMOVING REDUNDANT PREDICTED MICRORNAS AND FILTERING PREDICTED ONES BY PLANT CRITERIA
			#print 'perl '+dir+'/rm_redundant_meet_plant.pl '+output_dir+'/chromosome_length.txt '+infile+'_precursors.fa '+ infile+'_predictions '+ infile+'_nr_prediction '+infile+'_filter_P_prediction'
			os.system('nice -n 19 perl '+dir+'/rm_redundant_meet_plant.pl '+output_dir+'/chromosome_length.txt '+infile+'_precursors.fa '+ infile+'_predictions '+ infile+'_nr_prediction '+infile+'_filter_P_prediction') 			
		
			os.system ('mv '+output_dir+'/module_'+str(module_no)+'/miRNA_filter_P_prediction'+' '+output_dir+'/module_'+str(module_no)+'/miRNA_predictions')
		
			### making output that can be added to the MySQL table 	
			os.system('>'+output_dir+'/temp_folder/mysql_batch')
			infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
			os.system('nice -n 19 python '+script_dir+'/12af.py '+infile+' '+output_dir+'/module_'+str(module_no)+'/miRNA_predictions'+' '+output_dir+'/'+'check.txt '+' 2>'+output_dir+'/STDOUT.txt')
			os.system('nice -n 19 python '+script_dir+'/12t.py '+infile+' '+output_dir+'/'+'check.txt '+ 'miRNA_prediction_mirdeepP' +' >>'+output_dir+'/temp_folder/mysql_batch'+' 2>>'+output_dir+'/STDOUT.txt')
			os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/mysql_batch'+' 2>>'+output_dir+'/STDOUT.txt')
			os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/index.txt'+' 2>'+output_dir+'/STDOUT.txt')
	
			### move .txt file to module-4 folder
			os.system('mv '+output_dir+'/*.txt '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>>'+output_dir+'/STDOUT.txt')
				
			os.system('rm '+output_dir+'/module_'+str(module_no)+'/miRNA'+' 2>>'+output_dir+'/STDOUT.txt')	
			
	os.system('echo miRNA_predictions --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')
			
### funtion for module -4.2:	map the miRNA to mirBase
def miRNA_dataset_mapping(configuration,output_dir):
	import os
	t = time.time()
	os.system('echo miRBase_mapping --- starts >>time_stats.txt')
	module_no=4
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	### function for making a file redundant based on first column
	def make_redundant(infile):
		hash={}
		os.system('cp '+infile+' temp.txt')
		o=open(infile,'w')
		for line in open('temp.txt','r'):
			if line[0] in 'ATGCatgc':
				line=line.strip()
				token=line.split()
				hash[token[0]]=line
		for key in sorted(hash.iterkeys()):
			o.write(hash[key]+'\n')
		o.close()
		return infile

	if(len(configuration.mapping_infile)==0):
		infile=output_dir+'/module_'+str(module_no)+'/miRNA_predictions'
		os.system("awk 'BEGIN { FS = "+'"\t"'+"} {temp=$1; $1=$7; $7=temp} {print}' "+infile+ ' >'+infile+'_miRNA')
		infile=infile+'_miRNA'
	else:
		infile=configuration.mapping_infile	
	infile=make_redundant(infile)
	miRNA_dataset=configuration.miRNA_database
	#print 'nice -n 19 python '+script_dir+'/12d.py '+miRNA_dataset +' > '+miRNA_dataset+"_T"
	os.system('nice -n 19 python '+script_dir+'/12d.py '+miRNA_dataset +' > '+miRNA_dataset+"_T")
	#print 'nice -n 19 formatdb -i '+miRNA_dataset+"_T -p F -o T"
	os.system('nice -n 19 formatdb -i '+miRNA_dataset+"_T -p F -o T"+'\n')
	outfile=infile+'.fasta'
	#print 'nice -n 19 python '+script_dir+'/12c.py '+infile+' '+outfile
	os.system('nice -n 19 python '+script_dir+'/12c.py '+infile+' '+outfile)
	#print 'nice -n 19 blastall -p blastn -i '+outfile+'  -d '+miRNA_dataset+"_T -m 8 -o "+outfile+".blastout"
	os.system('nice -n 19 blastall -p blastn -i '+outfile+'  -d '+miRNA_dataset+"_T -m 8 -o "+outfile+".blastout")
	#print 'nice -n 19 python '+script_dir+'/12e.py '+ miRNA_dataset+"_T "+outfile+".blastout > "+outfile+".blastout.mismatch_count"
	os.system('nice -n 19 python '+script_dir+'/12e.py '+ miRNA_dataset+"_T "+outfile+".blastout > "+outfile+".blastout.mismatch_count")
	#print 'nice -n 19 python '+script_dir+'/12f.py '+infile+' '+ outfile+".blastout.mismatch_count >"+infile+".miRNA_mapped"
	os.system('nice -n 19 python '+script_dir+'/12f.py '+infile+' '+ outfile+".blastout.mismatch_count >"+infile+".miRNA_mapped")
	#os.system('head -n 2 '+infile+".miRNA_mapped")
	
	input=infile+".miRNA_mapped" # file should be tab seperated
	out=input+".miRNAfamily.txt"
	import linecache
	line = (linecache.getline(input, 1)).strip()
	column_no=str(len(line.split('\t')))
	os.system('>'+out)
	os.system('echo "Number of miRNA in database ">> '+out)
	count=int(commands.getoutput('echo `cat '+miRNA_dataset+'|wc -l`')) ### count of lines in fasta file/2
	os.system('echo '+ str(count/2) +' >>'+out)
	os.system('echo "Total number of miRNAs" >>'+out)
	os.system('nice -n 19 python '+script_dir+'/12a.py '+input +' '+column_no+' star|wc -l >>'+out)
	os.system('echo "Total number of miRNA families" >>'+out)
	os.system('nice -n 19 python '+script_dir+'/12a.py '+ input+' '+column_no+' star|sort -nr|uniq -c|wc -l>>'+out)
	os.system('echo "Total number of star sequences(also included in miRNA count)">>'+out)
	os.system('cat star >>'+out)
	
	os.system('echo "homologue for each family">>'+out)
	os.system('nice -n 19 python '+script_dir+'/12a.py '+ input+' '+ column_no+' star|sort -nr|uniq -c|sort -nr >test')
	os.system("awk '{ print $2"+'\t'+"$1 }' test >>"+out)

	os.system('rm '+output_dir+'/module_'+str(module_no)+'/miRNA_filter_P_prediction_miRNA'+' 2>'+output_dir+'/STDOUT.txt')
	
	os.system('echo miRBase mapping --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')

### funtion for module -5:	ta-siRNA prediction
def tasiRNA_prediction(configuration,output_dir):
	module_no=5
	import os
	t=time.time()
	os.system('echo ta-siRNA prediction --- >>time_stats.txt')
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	if(len(configuration.tasi_cluster_infile)==0):
		cluster_infile=output_dir+'/module_3'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+'.sorted.igv_cluster.bed'
	else:
		cluster_infile=configuration.tasi_cluster_infile
	
	if(len(configuration.tasi_mapping_infile)==0):
		mapping_infile=output_dir+'/module_2'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+'.sorted.igv'
	else:
		mapping_infile=configuration.tasi_mapping_infile
	os.system('nice -n 19 python '+script_dir+'/12y.py '+mapping_infile) 
	os.system('nice -n 19 perl '+script_dir+'/tasi.pl')
	genome_file=configuration.genome_file
	if(len(configuration.genome_file)==0):
		print "please provide correct reference file for ta-siRNA predictions"
	else:
		os.system('formatdb -i '+genome_file+' -p F -o T')
		os.system('nice -n 19 python '+script_dir+'/12z.py '+ cluster_infile+' '+genome_file+' |sort|uniq > '+output_dir+'/tasi.sh')
		os.system('nice -n 19 sh '+output_dir+'/tasi.sh >tasi_cluster.fa')

	os.system('rm sRNAmapping.txt'+' 2>'+output_dir+'/STDOUT.txt')
	os.system('rm tasi.sh'+' 2>'+output_dir+'/STDOUT.txt')
	#os.system('mv '+''+'*.sh '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
	os.system('mv '+''+'*.fa '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
	os.system('mv '+''+'*sRNA* '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')

	### making output that can be added to the MySQL table 	
	os.system('>'+output_dir+'/temp_folder/mysql_batch')
	infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
	os.system('nice -n 19 python '+script_dir+'/12ag.py '+infile+' '+output_dir+'/module_'+str(module_no)+'/p0.001_sRNA21nt_out.txt'+' '+output_dir+'/'+'check.txt '+' 2>'+output_dir+'/STDOUT.txt')
	os.system('nice -n 19 python '+script_dir+'/12t.py '+infile+' '+output_dir+'/'+'check.txt '+ 'ta-siRNAprediction' +' >>'+output_dir+'/temp_folder/mysql_batch'+' 2>>'+output_dir+'/STDOUT.txt')
	os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/mysql_batch'+' 2>>'+output_dir+'/STDOUT.txt')
	os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/index.txt'+' 2>>'+output_dir+'/STDOUT.txt')
	
	os.system('echo ta-siRNA prediction --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')
	
### funtion for module -6: adding annotation to mysql database
def add_annotation_to_mysql_database(configuration,output_dir):
	module_no=6
	import os
	t=time.time()
	os.system('echo annotation to mysql --- starts >>time_stats.txt')
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))

	if(len(configuration.infile)==0):
		infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
	else:
		infile=configuration.infile
	outfile=infile+'.fasta'
	os.system('nice -n 19 python '+script_dir+'/12c.py '+ infile +' '+ outfile+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt')
	for file in configuration.annotation_file:
		os.system('>'+output_dir+'/temp_folder/mysql_batch')
		### index the reference
		os.system("bowtie-build -q "+file+' '+file+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt')
		### map the reads
		os.system('nice -n 19 bowtie -q --threads '+str(configuration.no_of_cores)+' --seedmms '+str(configuration.mismatches)+' --maqerr '+str(configuration.max_maq_err)+' --seedlen '+str(configuration.seedlen)+' -k 100 --best --time --offbase 0 --sam '+ file +' -f '+ outfile+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt' +"| nice -n 19 perl -lane '"+"print $_ if ($F[3] > 0)' "+'|sort -nr -k 5,5 '+" > "+ infile+'.sam')
		### make file with sequences and annotations
		os.system('nice -n 19 python '+script_dir+'/12s.py '+infile+' '+ infile+'.sam '+ file + ' >'+output_dir+'/'+'out.txt'+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt')
		os.system('nice -n 19 python '+script_dir+'/12t.py '+infile+' '+output_dir+'/'+'out.txt '+ file +' >> '+output_dir+'/temp_folder/mysql_batch'+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt')
		os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/mysql_batch'+'\n')
		os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/index.txt'+'\n')
	
	if(configuration.add_genomic_region==True):
		infile=output_dir+'/module_2'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
		infile=infile+'.sorted.igv'
		os.system('nice -n 19 python '+script_dir+'/12h.py '+configuration.gtf_file+' '+infile +' |awk '+"'"+"{print $4 "+'"\t"'+" $NF}"+"'"+'>'+infile+'.region')
		exonic=0;intronic=0;overlap=0;intergenic=0;
		for line in open(infile+'.region','r'):
			line = line.strip()
			token=line.split()
			if (token[1]=='exonic'):
				exonic += 1
			if (token[1]=='intronic'):
				intronic += 1
			if (token[1]=='intergenic'):
				intergenic += 1
			if (token[1]== "exon-intron overlap"):
				overlap +=1
		sum=exonic+intronic+intergenic+overlap
		
		### making pichart for small RNA genomic distribution
		import matplotlib
		#matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		plt.pie([exonic/float(sum),intronic/float(sum),intergenic/float(sum),overlap/float(sum)],colors=((0.2023,1,0.8),(0.2023,0.8,1),(1,0.4745,0.3),(0.8509,0.4,1)),shadow=True,labels=['exonic','intronic','intergenic','overlapping'],pctdistance=0.6, labeldistance=1.1,autopct='%1.f%%',explode=(0.0,0.1,0.1,0.05))
		plt.title('Genomic Distribution of Small RNAs')
		plt.savefig(infile+'.region.png')
		
		o=open('sRNA_reads_Genomic_region_summary.txt','w')
		o.write( "number of exonic sequences: "+ str(exonic)+'\n')
		o.write( "number of intronic sequences: "+ str(intronic)+'\n')
		o.write( "number of intergenic sequences: "+ str(intergenic)+'\n')
		o.write( "number of exon-intron overlapping sequences: "+str(overlap)+'\n')
		o.close()
		
		
		os.system('>'+output_dir+'/temp_folder/mysql_batch')
		infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
		os.system('nice -n 19 python '+script_dir+'/12aa.py '+infile+' '+output_dir+'/module_2'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)+'.sorted.igv.region'+' '+output_dir+'/'+'check.txt ')
		os.system('nice -n 19 python '+script_dir+'/12t.py '+infile+' '+output_dir+'/'+'check.txt '+ 'genomic_region' +' >>'+output_dir+'/temp_folder/mysql_batch')
		os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/mysql_batch'+'\n')
		os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+output_dir+'/temp_folder/index.txt'+'\n')
		
		### move to module-6 folder
		os.system('mv '+output_dir+'/module_2/*region*'+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
		os.system('mv '+'*region*'+' '+str(date)+'_output'+'/module_'+str(module_no) + ' 2>'+output_dir+'/STDOUT.txt')
		
		os.system('echo annotation to mysql --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')

### funtion for module -7.1: making plot for expression values for libraries
def pattern_plot(configuration,output_dir):
	### generate a text file with 1's of size total_no_of_libraries*total_no_of_libraries
	import os
	t=time.time()
	os.system('echo expression plotting --- starts >>time_stats.txt')
	module_no=7
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)+'/plots'):
		os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)+'/plots')
	import os.path
	if not os.path.isfile(output_dir+'/module_1/'+'libraries_reference'):
		print 'library references are not provided'
	else:
		os.system('nice -n 19 python '+script_dir+'/12o.py '+str(len(configuration.lib))+' '+output_dir+'/module_1/'+'libraries_reference > '+output_dir+'/module_1/'+'compare_libraries')
	if(configuration.pattern_plot_file==''):
		for i in range(configuration.min_seq_size, configuration.max_seq_size+1,1):
			infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(i)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'
			os.system('nice -n 19 python '+script_dir+'/12w.py '+ infile+' '+str(configuration.score_regulation_cutoff)+' > '+infile+'_score_cutoff_'+str(configuration.score_regulation_cutoff))
			profile=infile+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
			#os.system('nice -n 19 python '+script_dir+'/spr.py '+profile+' >'+profile+'.spearman_values')
			os.system('nice -n 19 python '+script_dir+'/12p.py '+profile+' '+output_dir+'/module_1/'+'compare_libraries' +' '+output_dir+'/module_1/'+ 'libraries_reference '+profile+'_normalized.logScore '+profile+'_normalized.relativeFraction '+profile+'.spearman_values_log '+profile+'.spearman_values_norm ')
		infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'
		os.system('nice -n 19 python '+script_dir+'/12w.py '+ infile+' '+str(configuration.score_regulation_cutoff)+' > '+infile+'_score_cutoff_'+str(configuration.score_regulation_cutoff))
		profile=infile+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
		#print profile
		#os.system('nice -n 19 python '+script_dir+'/spr.py '+profile+' >'+profile+'.spearman_values')
		os.system('nice -n 19 python '+script_dir+'/12p.py '+profile+' '+output_dir+'/module_1/'+'compare_libraries' +' '+output_dir+'/module_1/'+ 'libraries_reference '+profile+'_normalized.logScore '+profile+'_normalized.relativeFraction '+profile+'.spearman_values_log '+profile+'.spearman_values_norm ')
	else:
		profile=configuration.pattern_plot_file
		#print profile
		os.system('nice -n 19 python '+script_dir+'/12p.py '+profile+' '+output_dir+'/module_1/'+'compare_libraries' +' '+output_dir+'/module_1/'+ 'libraries_reference '+profile+'_normalized.logScore '+profile+'_normalized.relativeFraction '+profile+'.spearman_values_log '+profile+'.spearman_values_norm ')
		#os.system('nice -n 19 python '+script_dir+'/spr.py '+profile+' >'+profile+'.spearman_values')
	os.system('mv '+output_dir+'/module_1/profiles'+'/*png '+output_dir+'/module_'+str(module_no)+'/plots/')
	os.system('mv '+output_dir+'/module_1/profiles'+'/*.spearman_values* '+output_dir+'/module_'+str(module_no))
	os.system('rm '+output_dir+'/module_1/profiles'+'/*.relativeFraction* ')
	os.system('rm '+output_dir+'/module_1/profiles'+'/*.logScore* ')
	os.system('echo expression plotting --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')
	
### funtion for module -1.4: making plots for size distribution
def size_distribution(configuration,output_dir):
	import os
	module_no=1
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)+'/plots'):
		os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)+'/plots')

	if(len(configuration.infile)==0):
		infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
	else:
		infile=configuration.infile
	#print "number of sequence: ",int(commands.getoutput('echo `cat '+infile+'|wc -l`'))-1
	os.system('nice -n 19 python '+script_dir+'/12j.py '+ infile +' > '+date+'_size_distribution_'+infile.split('/')[-1])
	os.system('mv '+output_dir+'/module_1/profiles'+'/*png '+output_dir+'/module_'+str(module_no))
	os.system('mv '+date+'_size_distribution_'+infile.split('/')[-1]+' '+output_dir+'/module_'+str(module_no)+'/plots/')
	
	
### funtion for module -8.1: find pattern in first three nucleotides, morkov probabilities
def find_pattern(configuration,output_dir):
	import os
	t=time.time()
	os.system('echo find pattern --- starts >>time_stats.txt')
	module_no=8
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	if(len(configuration.infile_for_pattern)==0):
		infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
	else:
		infile=configuration.infile_for_pattern
	if(len(configuration.ref_file_for_pattern)==0):
		print 'please provide reference genome file'
	else:
		### calculate probabability of A,T,G,C
		nuc='A','C','G','T';count_nuc={'A':0,'C':0,'G':0,'T':0}
		for line in open(configuration.ref_file_for_pattern,'r'):
			line=line.strip()
			for i in range(4):
				count_nuc[nuc[i]] += line[0].count(nuc[i])
		o=open('prob','w')
		for n in nuc:
			o.write(str(n)+'\t'+str(count_nuc[n]/float(sum(count_nuc.values())))+'\n')
		o.close()
		os.system('nice -n 19 perl '+script_dir+'/MarkovModel.pl '+ infile +' '+'>'+output_dir+'/markov_prob.txt')
		os.system('mv '+output_dir+'/markov_prob.txt'+' '+output_dir+'/module_'+str(module_no))
		os.system('echo find pattern --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')
	
### funtion for module -8.2:	map the sequences with the pattern
def map_pattern(configuration,output_dir):
	import os
	t=time.time()
	os.system('echo map pattern --- starts >>time_stats.txt')
	module_no=8
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	if(len(configuration.infile_for_pattern_map)==0):
		infile=output_dir+'/module_1/profiles'+'/'+date+'_profile_'+str(configuration.min_seq_size)+'_'+str(configuration.max_seq_size)+'_sorted'+'.cut-'+str(configuration.cutoff)+'_normalized'+'_score_cutoff_'+str(configuration.score_regulation_cutoff)
	else:
		infile=configuration.infile_for_pattern_map
	outfile=infile+'.fasta'
	file=configuration.genome_file
	os.system('nice -n 19 python '+script_dir+'/12c.py '+ infile +' '+ outfile)
	os.system("bowtie-build -q "+file+' '+file+' '+' 2>>'+output_dir+'/module_'+str(module_no)+'/pattern_map.log')
	os.system('nice -n 19 bowtie -q --threads '+str(configuration.no_of_cores)+' --seedmms '+str(configuration.mismatches)+' --maqerr '+str(configuration.max_maq_err)+' --seedlen '+str(configuration.seedlen)+' -k 30 --best --time --offbase 0 --sam '+ file +' -f '+ outfile +' 2>>'+output_dir+'/module_'+str(module_no)+'/pattern_map.log'+"| nice -n 19 perl -lane '"+"print $_ if ($F[3] > 0)' > "+ outfile+'.sam'+' 2>'+output_dir+'/temp_folder/ERROR_OUT.txt')
	x=int(commands.getoutput('echo `cat '+outfile+'.sam'+'|wc -l`'))
	y=int(commands.getoutput('echo `cut -f 1 '+outfile+'.sam'+'|sort| uniq -c| wc -l`'))
	z=int(commands.getoutput('echo `cat '+outfile+'|wc -l`'))/2
	
	o=open(output_dir+'/module_'+str(module_no)+'/pattern_read_stats.txt','w')
	o.write("Total number of unique sequences: "+str(z)+'\n')
	o.write("Total number of places mapped on reference: "+str(x)+'\n')
	o.write("Total number of unique sequences mapped: "+str(y)+'\n')
	o.write("Percentage of mapped unique sequences: "+str(round(y*100/(z),2))+'%\n')
	
	if(len(configuration.pattern)!=0):
		infile=outfile
		os.system("nice -n 19 python "+ script_dir+"/12b.py "+infile+" "+configuration.pattern+' > '+infile+'_'+configuration.pattern)
		outfile=infile+'_'+configuration.pattern
		os.system("bowtie-build -q "+file+' '+file+' '+' 2>>'+output_dir+'/module_'+str(module_no)+'/pattern_map.log')
		os.system('nice -n 19 bowtie -q --threads '+str(configuration.no_of_cores)+' --seedmms '+str(configuration.mismatches)+' --maqerr '+str(configuration.max_maq_err)+' --seedlen '+str(configuration.seedlen)+' -k 30 --best --time --offbase 0 --sam '+ file +' -f '+ outfile +' 2>>'+output_dir+'/module_'+str(module_no)+'/pattern_map.log' +"| nice -n 19 perl -lane '"+"print $_ if ($F[3] > 0)' > "+ outfile+'.sam'+' 2>>'+output_dir+'/temp_folder/ERROR_OUT.txt')
		x=int(commands.getoutput('echo `cat '+outfile+'.sam'+'|wc -l`'))
		y=int(commands.getoutput('echo `cut -f 1 '+outfile+'.sam'+'|sort| uniq -c| wc -l`'))
		z=int(commands.getoutput('echo `cat '+outfile+'|wc -l`'))/2
		o.write("total number of unique sequences with pattern: "+str(z)+'\n')
		o.write("total number of places mapped on reference with pattern: "+str(x)+'\n')
		o.write("total number of unique sequences mapped with pattern: "+ str(y)+'\n')
		o.write("Percentage of mapped unique sequences : "+str(round(y*100/(z),2))+'%\n')
		o.close()
		
	### move to module -8
	os.system('mv '+output_dir+'/module_1/profiles'+'/*GAA '+output_dir+'/module_'+str(module_no))
	os.system('rm '+output_dir+'/module_1/profiles'+'/*GAA* ')
	os.system('rm '+output_dir+'/module_1/profiles'+'/*.sam ')
	
	os.system('echo map_pattern --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')

### funtion for module -9:	make mySQL queries
def mySQL_query(configuration,output_dir):
	import os
	t=time.time()
	os.system('echo mysql query --- starts >>time_stats.txt')
	module_no=9
	if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no)):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/module_'+str(module_no))
	if(len(configuration.query_infile)>0):
		os.system('mysql -h '+configuration.server+' -u '+configuration.user+' -p '+ configuration.database+' --password='+str(configuration.password)+' <'+configuration.query_infile+' >mysql_query_output')
	
	### move to module-9
	os.system('mv '+'mysql_query_output'+' '+output_dir+'/module_'+str(module_no))
	os.system('echo mysql query --- ends time taken: '+str(time.time()-t)+' >>time_stats.txt')
	
if __name__ == "__main__":


	if (configure_run=="N" or configure_run=="n"):
		print "your output will be in folder name: "+str(date)+'_output'
		### makig a output folder
		if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output')
		if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'+'/temp_folder'):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/temp_folder')
		output_dir=os.getcwd()+'/'+str(date)+'_output'
		os.system('>'+output_dir+'/temp')
		
		
				
		'''Module-1.1'''
		print "%-40s%-10s" %("1.1  :Raw data filtering",configuration_12.Raw_data_filtering)

		r=configuration_12.Raw_data_filtering
		### check if filtering needs to be done
		if (r==False) or (r=='False') or (r=='F') or (r=='f'):
			r=configuration_12.Make_profile_files	
			if(r==True) or (r=='True') or (r=='T') or (r=='t'):
				path=output_dir+'/module_1'
				filtered=False
				for infile in glob.glob( os.path.join(path, '*')):
					if(re.search("repeat", infile)):
						filtered=True
				if(filtered!=True):
					configuration_12.Raw_data_filtering=True
					configuration_12.repeat_database='empty_file.txt'
	
		
				
		r=configuration_12.Raw_data_filtering
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			Raw_data_filtering(configuration_12,output_dir) 
		
		'''Module-1.2'''
		print "%-40s%-10s" %('1.2  :Make_profile_files:',configuration_12.Make_profile_files)
		r=configuration_12.Make_profile_files
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			Make_profile_files(configuration_12,output_dir)
		
		'''Module-1.3'''
		print "%-40s%-10s" %("1.3  :Add_output_to_mysql_database:",configuration_12.Add_output_to_mysql_database)
		r=configuration_12.Add_output_to_mysql_database
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			add_output_to_mysql_database(configuration_12,output_dir)
				
		'''Module-2'''
		print "%-40s%-10s" %("2    :Genome_mapping:",configuration_12.Genome_mapping)
		r=configuration_12.Genome_mapping
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			genome_mapping(configuration_12,output_dir)			
	
		
		'''Module-3.1'''
		print "%-40s%-10s" %("3.1  :Cluster_prediction:",configuration_12.Cluster_prediction)
		r=configuration_12.Cluster_prediction
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			cluster_prediction(configuration_12,output_dir)
			
		'''Module-3.2'''
		'''
		print "%-40s%-10s" %("3.x  :Calculate_unique_sequences_in_cluster:",configuration_12.Calculate_unique_sequences_in_cluster)
		r=configuration_12.Calculate_unique_sequences_in_cluster
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			calculate_unique_sequences_in_cluster(configuration_12,output_dir)
			
		"""Module-3.3"""
		print "%-40s%-10s" %("3.x  :Regulated_sequences_in_cluster:",configuration_12.Regulated_sequences_in_cluster)
		r=configuration_12.Regulated_sequences_in_cluster
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			regulated_sequences_in_cluster(configuration_12,output_dir)
		'''
			
		'''Module-3.2'''
		print "%-40s%-10s" %("3.2  :Cluster_genomic_region_prediction:",configuration_12.Cluster_genomic_region_prediction)
		r=configuration_12.Cluster_genomic_region_prediction
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			cluster_genomic_region_prediction(configuration_12,output_dir)
						
		'''Module-4.1'''
		print "%-40s%-10s" %("4.1  :MiRNA_prediction:",configuration_12.MiRNA_prediction)
		r=configuration_12.MiRNA_prediction
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			miRNA_prediction(configuration_12,output_dir)
		
		'''Module-4.2'''
		print "%-40s%-10s" %("4.2  :MiRNA_dataset_mapping:",configuration_12.MiRNA_dataset_mapping)
		r=configuration_12.MiRNA_dataset_mapping
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			miRNA_dataset_mapping(configuration_12,output_dir)
					
		'''Module-5'''
		print "%-40s%-10s" %("5    :ta-siRNA prediction:",configuration_12.TasiRNA_prediction)
		r=configuration_12.TasiRNA_prediction
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			tasiRNA_prediction(configuration_12,output_dir)
			
			
		'''Module-6'''
		print "%-40s%-10s" %("6    :Add_annotation_to_mysql_database:",configuration_12.Add_annotation_to_mysql_database)
		r=configuration_12.Add_annotation_to_mysql_database
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			add_annotation_to_mysql_database(configuration_12,output_dir)

		'''Module-7.1'''
		print "%-40s%-10s" %("7    :Pattern_plotting:",configuration_12.Pattern_plot)
		r=configuration_12.Pattern_plot
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			pattern_plot(configuration_12,output_dir)
		
				
		'''Module-8.1'''
		print "%-40s%-10s" %("8.1  :Find_pattern:",configuration_12.Find_pattern)
		r=configuration_12.Find_pattern
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			find_pattern(configuration_12,output_dir)
		
		'''Module-8.2'''
		print "%-40s%-10s" %("8.2  :Map_pattern:",configuration_12.Map_pattern)
		r=configuration_12.Map_pattern
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			map_pattern(configuration_12,output_dir)
		
		'''Module-9'''
		print "%-40s%-10s" %("9    :MySQL_query:",configuration_12.MySQL_query)
		r=configuration_12.MySQL_query
		if(r==True) or (r=='True') or (r=='T') or (r=='t'):
			mySQL_query(configuration_12,output_dir)

			
		
		
	else:	
		### makig a output folder
		print "your output will be in folder name: "+str(date)+'_output'
		if not os.path.exists(os.getcwd()+'/'+str(date)+'_output'):
			os.mkdir(os.getcwd()+'/'+str(date)+'_output')
			os.mkdir(os.getcwd()+'/'+str(date)+'_output'+'/temp_folder')
		output_dir=os.getcwd()+'/'+str(date)+'_output'
		os.system('>'+output_dir+'/temp')
		
		'''Module-1.1'''
		print "%-40s%-10s" %("1.1  :Raw data filtering:",configuration_12i.Raw_data_filtering)

			
		r=configuration_12i.Raw_data_filtering
		### check if filtering needs to be done
		if (r==False) or (r=='False') or (r=='F') or (r=='f'):
			r=configuration_12i.Make_profile_files	
			if(r==True) or (r=='True') or (r=='T') or (r=='t'):
				path=output_dir+'/module_1'
				filtered=False
				for infile in glob.glob( os.path.join(path, '*')):
					if(re.search("repeat_filtered", infile)):
						filtered=True
				if(filtered!=True):
					configuration_12i.Raw_data_filtering=True
					configuration_12i.repeat_database='empty_file.txt'
	
		
		
		if(configuration_12i.Raw_data_filtering==True):
			Raw_data_filtering(configuration_12i,output_dir) 
		
		'''Module-1.2'''
		print "%-40s%-10s" %('1.2  :Make_profile_files:',configuration_12i.Make_profile_files)
		
		if(configuration_12i.Make_profile_files==True):
			Make_profile_files(configuration_12i,output_dir)
		
		'''Module-1.3'''
		print "%-40s%-10s" %("1.3  :Add_output_to_mysql_database:",configuration_12i.Add_output_to_mysql_database)	
		if(configuration_12i.Add_output_to_mysql_database==True):
			add_output_to_mysql_database(configuration_12i,output_dir)
				
		'''Module-2'''
		print "%-40s%-10s" %("2    :Genome_mapping:",configuration_12i.Genome_mapping)
		if(configuration_12i.Genome_mapping==True):
			genome_mapping(configuration_12i,output_dir)			
	
		
		'''Module-3.1'''
		print "%-40s%-10s" %("3.1  :Cluster_prediction:",configuration_12i.Cluster_prediction)
		if(configuration_12i.Cluster_prediction==True):
			cluster_prediction(configuration_12i,output_dir)
			
		'''
		"""Module-3.2"""
		print "%-40s%-10s" %("3.x  :Calculate_unique_sequences_in_cluster:",configuration_12i.Calculate_unique_sequences_in_cluster)
		if(configuration_12i.Calculate_unique_sequences_in_cluster==True):
			calculate_unique_sequences_in_cluster(configuration_12i,output_dir)
			
		"""Module-3.3"""
		print "%-40s%-10s" %("3.x  :Regulated_sequences_in_cluster:",configuration_12i.Regulated_sequences_in_cluster)
		if(configuration_12i.Regulated_sequences_in_cluster==True):
			regulated_sequences_in_cluster(configuration_12i,output_dir)
		'''
			
		'''Module-3.4'''
		print "%-40s%-10s" %("3.2  :Cluster_genomic_region_prediction:",configuration_12i.Cluster_genomic_region_prediction)
		if(configuration_12i.Cluster_genomic_region_prediction==True):
			cluster_genomic_region_prediction(configuration_12i,output_dir)
						
		'''Module-4.1'''
		print "%-40s%-10s" %("4.1  :MiRNA_prediction:",configuration_12i.MiRNA_prediction)
		if(configuration_12i.MiRNA_prediction==True):
			miRNA_prediction(configuration_12i,output_dir)
		
		'''Module-4.2'''
		print "%-40s%-10s" %("4.2  :MiRNA_dataset_mapping:",configuration_12i.MiRNA_dataset_mapping)
		if(configuration_12i.MiRNA_dataset_mapping==True):
			miRNA_dataset_mapping(configuration_12i,output_dir)
					
		'''Module-5'''
		print "%-40s%-10s" %("5    :ta-siRNA prediction:",configuration_12i.TasiRNA_prediction)
		if(configuration_12i.TasiRNA_prediction==True):
			tasiRNA_prediction(configuration_12i,output_dir)
			
					
		'''Module-6'''
		print "%-40s%-10s" %("6    :Add_annotation_to_mysql_database:",configuration_12i.Add_annotation_to_mysql_database)
		if(configuration_12i.Add_annotation_to_mysql_database==True):
			add_annotation_to_mysql_database(configuration_12i,output_dir)

		'''Module-7'''
		print "%-40s%-10s" %("7    :Pattern_plotting:",configuration_12i.Pattern_plot)
		if(configuration_12i.Pattern_plot==True):
			pattern_plot(configuration_12i,output_dir)
		
				
		'''Module-8.1'''
		print "%-40s%-10s" %("8.1  :Find_pattern:",configuration_12i.Find_pattern)
		if(configuration_12i.Find_pattern==True):
			find_pattern(configuration_12i,output_dir)
		
		'''Module-8.2'''
		print "%-40s%-10s" %("8.2  :Map_pattern:",configuration_12i.Map_pattern)
		if(configuration_12i.Map_pattern==True):
			map_pattern(configuration_12i,output_dir)
		
		'''Module-9'''
		print "%-40s%-10s" %("9    :MySQL_query:",configuration_12i.MySQL_query)
		if(configuration_12i.MySQL_query==True):
			mySQL_query(configuration_12i,output_dir)
