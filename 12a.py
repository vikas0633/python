### script for making summary on miRNA hits
import sys,re

column=int(sys.argv[2])
star=open(sys.argv[3],'w')
st=0
### open input file
for line in open(sys.argv[1],'r'):
	line=line.strip()
	token=line.split('\t')
	sta=token[column-1].split('*')
	miRNA=token[column-1].replace('.','-')
	if(len(sta)>1):
		if(len(sta[1])==0):
			st += 1 ###count number of stars
	tokens=miRNA.split('-')
	try:
		mir=re.sub("[a-zA-Z]",'',tokens[1])
		mir=mir.replace('*','')
		if(len(mir)>0):
			print("MIR"+str(int(mir)))
	except:
		continue
		
star.write(str(st))
star.write('\n')
star.close()