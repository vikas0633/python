### shorting position on chromosome and writng to a new file all inormation we need
### Uses igv files as input
import time,numpy
import sys, os
#sys.setrecursionlimit(5000)

work_dir=os.getcwd()
lib_path = os.path.abspath(os.getcwd())
script_dir=os.getcwd()+'/'+'scripts'
sys.path.append(lib_path)
from configuration_12 import genotype
genotype=list(genotype)

        
def find_cluster(abud_cutoff,fi,fo,gap,mur):
        f = open(fi,'r')
        out = open(fo,'w')
        out.write(str("Chrom"))
        out.write("\t")
        out.write(str("ChromStart"))
        out.write("\t")
        out.write(str("ChromEnd"))
        out.write("\t")
        out.write(str("Name"))
        out.write("\t")
        out.write(str("Abundance"))
        out.write("\n")
        out.close()
        ### count different size sRNA
        size = {}
        i= -1 
        sRNA = {}
        previous = -10000
        s = 0
        abundance = 0
        cluster = {}
        cluster_sRNA = {}
        in_cluster = False
        ## reading lines
        for line in f.readlines():
                i += 1
                line = line.strip()
                token = line.split('\t')
                try:
                        genome_position = float(token[1])
                except:
                        continue
                sRNA[i] = token[3]
                chro = token[0]
                ab = 0
                for j in range(4,4+len(genotype)):
                    ab += float(token[j])
                if(i >=1):
                        #prfloat(genome_position,previous)
                        if((float(genome_position) - float(previous)) <= gap):
                                if(in_cluster == False):
                                        start = previous
                                        name = nam
                                        cluster_sRNA[0] = 1
                                        abundance = abund
                                        in_cluster = True
                                token=sRNA[i].split('_')
                                abundance += ab
                                #prfloat(abundance,token[3],sRNA[position])
                                cluster_sRNA[0] += 1
                                name += ','
                                name += sRNA[i]
                        if((float(genome_position) - (float(previous)) > gap)&(in_cluster==True)):
                                name = str(cluster_sRNA[0])+','+name
                                in_cluster = False
                                if(((cluster_sRNA[0])>=mur)&( abundance >=abud_cutoff)):
                                        end = previous+25
                                        out = open(fo,'a')
                                        out.write(str(chro))
                                        out.write('\t')
                                        out.write(str(int(start)))
                                        out.write('\t')
                                        out.write(str(int(end)))
                                        out.write('\t')
                                        out.write(str(name))
                                        out.write('\t')
                                        out.write(str(int(abundance)))
                                        out.write('\n')
                                        out.close()
                                        names = name.split(',')
                                        for t in range(1,len(names)-1):
                                        	if len(names[t]) in size:
                                        		size[len(names[t])] += 1
                                        	else:
                                        		size[len(names[t])] = 0
                                        #prfloat(chro,start,end,name,abundance)
                                        s += 1
                previous = genome_position
                nam = sRNA[i]
                token=sRNA[i].split('_')
                abund = ab
        print "abundance cut-off: "+str(abud_cutoff)
        print "max gap length: "+str(gap)
        print "min number fo unique reads: "+str(mur)
        print "number of clusters: "+str(s)
        print "size fractionated read counts in clusters:",size
        
def find_cluster_size(abud_cutoff,fi,fo,seq_size,gap,mur):
	f = open(fi,'r')
	out = open(str(fo)+'_'+str(seq_size)+".txt",'w')
	out.write(str("Chrom"))
	out.write("\t")
	out.write(str("ChromStart"))
	out.write("\t")
	out.write(str("ChromEnd"))
	out.write("\t")
	out.write(str("Name"))
	out.write("\t")
	out.write(str("Abundance"))
	out.write("\n")
	out.close()
	### count different size sRNA
	size = {}
	i= 0 
	sRNA = {}
	previous = -10000
	s = 0
	abundance = 0
	cluster = {}
	cluster_sRNA = {}
	in_cluster = False
	## reading lines
	for line in f.readlines():
		line = line.strip()
		token = line.split('\t')
		try:
				genome_position = float(token[1])
		except:
				continue
		if(len(token[3])==int(seq_size)): ###only read with specific size
			sRNA[i] = token[3]
			chro = token[0]
			ab = 0
			for j in range(4,4+len(genotype)):
					ab += float(token[j])
			if(i >=1):
					#prfloat(genome_position,previous)
					if((float(genome_position) - float(previous)) <= gap):
							if(in_cluster == False):
									start = previous
									name = nam
									cluster_sRNA[0] = 1
									abundance = abund
									in_cluster = True
							token=sRNA[i].split('_')
							abundance += ab
							#prfloat(abundance,token[3],sRNA[position])
							cluster_sRNA[0] += 1
							name += ','
							name += sRNA[i]
					if((float(genome_position) - (float(previous)) > gap)&(in_cluster==True)):
							name = str(cluster_sRNA[0])+','+name
							in_cluster = False
							if(((cluster_sRNA[0])>=mur)&( abundance >=abud_cutoff)):
									end = previous+25
									out = open(str(fo)+'_'+str(seq_size)+".txt",'a') ### difference size different file ??
									out.write(str(chro))
									out.write('\t')
									out.write(str(int(start)))
									out.write('\t')
									out.write(str(int(end)))
									out.write('\t')
									out.write(str(name))
									out.write('\t')
									out.write(str(int(abundance)))
									out.write('\n')
									out.close()
									names = name.split(',')
									for t in range(1,len(names)-1):
										if len(names[t]) in size:
											size[len(names[t])] += 1
										else:
											size[len(names[t])] = 0

									#prfloat(chro,start,end,name,abundance)
									s += 1
			previous = genome_position
			nam = sRNA[i]
			token=sRNA[i].split('_')
			abund = ab
			i += 1
	print "abundance cut-off: "+str(abud_cutoff)
	print "max gap length: "+str(gap)
	print "min number fo unique reads: "+str(mur)
	print("number of clusters: "+str(s))
	print "size:number_of_reads =", size 
""" for plotting size distribution         
	import numpy
	y = [0]*6
	x = [0]*6
	z = 0
	for element in size:
		x[z] = element -0.5
		y[z] = size[element]
		z += 1
	from pylab import *
	fig = figure()
	ax = fig.add_subplot(111,polar = False)
	ax.bar(x, y)
	show()
"""

def file_call(fi,fo,abu,gap,mur):
	start = time.clock()
	abud_cutoff = abu
	find_cluster(abud_cutoff,fi,fo,gap,mur)
	lap = time.clock() - start
        
def file_call_size(fi,fo,seq_size,abu,gap,mur):
	start = time.clock()
	abud_cutoff = abu
	find_cluster_size(abud_cutoff,fi,fo,int(seq_size),gap,mur)
	lap = time.clock() - start

        
### reading commond line

fi = sys.argv[1]
fo = sys.argv[2]
abu=float(sys.argv[4])
### call function to calculate clusters

### define cluster based on different sequence size
divide_by_size = sys.argv[3]

### max gap between two reads
gap=int(sys.argv[7])
min_unique_read=int(sys.argv[8])
mur=min_unique_read
print "size fractionation: ",(divide_by_size)

if (divide_by_size == 'True'):
	for seq_size in range(int(sys.argv[5]),int(sys.argv[6])+1):
		file_call_size(fi,fo,seq_size,abu,gap,mur)
	file_call(fi,fo,abu,gap,mur)
else:
	file_call(fi,fo,abu,gap,mur)
