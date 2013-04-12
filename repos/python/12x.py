###12x.py - script for making chromosome length file for miRNA predictions
import sys 




def cal_len(infile):
	### make hash for cds file
	cds={}
	key=''
	seq=''
	for line in open(infile,'r'):
		line=line.strip()
		token=line.split()
		if(line[0]=='>'):
			if(len(key)!=0):
				print key+'\t'+str(len(seq))
				seq=''
			key=token[0][1:]
		else:
			seq+=line
	print key+'\t'+str(len(seq))

if __name__ == "__main__":
	cal_len(sys.argv[1])