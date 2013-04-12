### script for making table the library wise abundances
import sys, os, glob

work_dir=os.getcwd()
lib_path = os.path.abspath(os.getcwd())
script_dir=os.getcwd()+'/'+'scripts'
sys.path.append(lib_path)
from configuration_12 import genotype



string=sys.argv[2]
store=[]
for line in open(sys.argv[1],'r'):
	line=line.strip()
	token=line.split(' ')
	store.append(token[1])
	
for i in range(1, len(genotype)+1):
	string += ('\t'+store[i])
	
print string