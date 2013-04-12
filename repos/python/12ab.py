### script for making header for abundance table
import os,sys,glob

work_dir=os.getcwd()
lib_path = os.path.abspath(os.getcwd())
script_dir=os.getcwd()+'/'+'scripts'
sys.path.append(lib_path)
from configuration_12 import genotype
genotype=list(genotype)


string='Size'

for i in genotype:
	string += ('\t'+i)
	
print string