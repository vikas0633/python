#12r.py - /plant/2011_week43/scripts -  for making sql batch script
import sys

file_name=sys.argv[1]
dir=sys.argv[2]
file=open(sys.argv[1],'r')
out=open(sys.argv[3],'w')

### making Table
first_line=False
for line in open(sys.argv[1],'r'):
	line = line.strip()
	if (len(line)>0):
		if(first_line==False):
			token=line.split()
			first_line=True
			string="`"+str(token[0])+"`"+" VARCHAR(100), "
			header="`"+str(token[0])+"`"+", "
			for i in range(1,len(token)-1):
				string += "`"+str(token[i])+"`"+" FLOAT "+","
			string += "`"+str(token[len(token)-1])+"`"+" FLOAT "
			for j in range(1,len(token)-1):
				if(j < 14):
					header += "`"+str(token[j])+"`"+", "
				else:
					continue
			header += "`"+str(token[len(token)-1])+"`"
		
file_name=(file_name.split('/'))[-1]
file_name=file_name.replace('_len','')
file_name=file_name.replace('normalized_score_cutoff','score')
string = str("CREATE TABLE "+'`'+file_name+'` '+ "("+string+");")
print(string)
print str("DROP TABLE IF EXISTS temp;")
print str("DROP TABLE IF EXISTS test;")


### loading data
string = "LOAD DATA LOCAL INFILE "+"'"+dir+'/'+"check.txt"+"'"+" INTO TABLE "+"`"+file_name+"`"+";"

print(string)

### creating index
string = str("CREATE INDEX "+'`'+file_name+"_index"+'` '+ " ON "+"`"+file_name+"`"+"("+header+");")
out.write(str(string+'\n'))
print(string)