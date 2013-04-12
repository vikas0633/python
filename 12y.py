## to make a compatible input for chen's algorithm
### from igv files
import sys
output = open('sRNAmapping.txt','w')
f = open(sys.argv[1],'r')
for line in f.readlines():
        line = line.strip()
        token = line.split()
        try:
                if(len(token[3])==21):
                       output.write(str(token[0]+'\t'+token[1]+'\t'+'1'+'\t'+token[3]+'\n'))
        except:
                print("passing")
f.close()
output.close()
