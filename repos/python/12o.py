#12o.py - /plant/2011_week42/20111019 - ### generate a text file with 0's of size total_no_of_libraries*total_no_of_libraries


import os,sys


work_dir=os.getcwd()
lib_path = os.path.abspath(os.getcwd())
script_dir=os.getcwd()+'/'+'scripts'
sys.path.append(lib_path)
from configuration_12 import genotype
genotype=list(genotype)

lib=int(sys.argv[1])

### print header
line='#\t'
for i in range(0,lib-1):
	line += genotype[i]+"\t"
line += genotype[i+1]+'\t'
line += 'library_sum'+'\t'
line += 'score'+'\n'
print (line)
for i in range(0,lib):
	line= genotype[i]+"\t"
	for j in range(0,i):
		line += (str(1)+'\t')
	line += (str(1)+'\n')
	print line
line= 'library_sum'+"\t"
for j in range(0,lib):
	line += (str(1)+'\t')
line += (str(1)+'\n')
print line
line= 'score'+"\t"
for j in range(0,lib):
	line += (str(1)+'\t')
line += (str(1)+'\t')
line += (str(1)+'\n')
print line

### open file for writing reference libraries
o=open(sys.argv[2],'w')

line='#'+'library'+'\t'+'reference'+'\n'
o.write(line)
for j in range(0,lib):
	if(j%2==0):
		o.write(str(j+1)+'\t'+str(0)+'\n')
	else:
		o.write(str(j+1)+'\t'+str(0)+'\n')
o.write(str('sum')+'\t'+str(0)+'\n')
o.write(str('score')+'\t'+str(0)+'\n')
