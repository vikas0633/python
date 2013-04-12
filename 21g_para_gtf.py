###21g_para_gtf.py - script for calculating exon/intron/transcripts lengths - /Users/vikas0633/Desktop/script/python 
### /Users/vikas0633/Desktop/

import sys, getopt,os, re
import numpy as np

### get the options 
def options(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:g:",["gtf=",'genome='])
	except getopt.GetoptError:
		print ''' 
	nice -n 19 python 21g_para_gtf.py 
    		-i <gtf file name> 
    		-g <genome file, optional>
    		'''
    		sys.exit(2)
	### default values
	gfile='';
	folder=os.getcwd()
	for opt, arg in opts:
		if opt == '-h':
			print ''' 
	nice -n 19 python 21g_para_gtf.py 
    		-i <gtf file name> 
    		-g <genome file, optional>
    		'''
			sys.exit()
		elif opt in ("-i", "--gtf"):
			infile = arg
		elif opt in ("-g", "--genome"):
			gfile = arg			
	return infile, gfile


### searching keywords in the memory

def search_file(word):
	word='\t'+word+'\t'
	index = re.search(word,text)
	if index :
		return True
	else:
		return False
	

### process data

def process_data(file):
	exon=[]
	intron=[]
	exon_all=[]
	intron_all=[]
	chr=[]
	genes=[]
	genes_all=[]
	last_gene=''
	for line in open(file,'r'):
		line = line.strip()
		if len(line)>0:
			token = line.split('\t')
			if len(token) > 4:
				### check if line is exon
				if token[2] == "exon":
					gene=line.split('"')[1]
					transcript = line.split('"')[3]
					if token[0] not in chr:
						exon.append([])
						intron.append([])
						genes.append([])
						chr.append(token[0])
						exon[len(exon)-1].append(abs(int(token[3])-int(token[4])))
						exon_all.append(abs(int(token[3])-int(token[4])))
					else:
						exon[len(exon)-1].append(abs(int(token[3])-int(token[4])))
						exon_all.append(abs(int(token[3])-int(token[4])))
						if (last_transcript == transcript) & (last_chr==token[0]):
							intron[len(intron)-1].append(abs(int(token[3])-last_end))
							intron_all.append(abs(int(token[3])-last_end))
					### check for gene length
					if last_gene != gene:
						if last_gene == '':
							start = int(token[3])
						else:
							genes[len(genes)-1].append(abs(long_end-start))
							genes_all.append(abs(long_end-start))
						start = int(token[3])
						long_end=0
					last_gene=gene
					last_chr=token[0]
					last_end=int(token[4])
					last_transcript=transcript
					### to get longest coverage of gene by all possible alternative splicing
					if long_end < last_end:
						long_end = last_end
	genes[len(genes)-1].append(abs(long_end-start))
	genes_all.append(abs(long_end-start))	
	print "Features       %20s %30s %30s %30s %30s %30s" %('Count','Total length','Average length','longest','shortest','avg per gene')
#	print "Gene           %20s %30s %30s %30s %30s %30s" %(str(len(genes_all)),str(np.sum(genes_all)),str(round(np.average(genes_all),2)),str(np.max(genes_all)),str(np.min(genes_all)),str(round(len(genes_all)/float(len(genes_all)),2)))
	print "Exon           %20s %30s %30s %30s %30s %30s" %(str(len(exon_all)),str(np.sum(exon_all)),str(round(np.average(exon_all),2)),str(np.max(exon_all)),str(np.min(exon_all)),str(round(len(exon_all)/float(len(genes_all)),2)))
	print "Intron         %20s %30s %30s %30s %30s %30s" %(str(len(intron_all)),str(np.sum(intron_all)),str(round(np.average(intron_all),2)),str(np.max(intron_all)),str(np.min(intron_all)),str(round(len(intron_all)/float(len(genes_all)),2)))
	
	count_genes=[]
	for i in genes:
		count_genes.append(len(i)) 
	
	count_exon=[]
	for i in exon:
		count_exon.append(len(i))
		
	count_intron=[]
	for i in intron:
		count_intron.append(len(i))
	
	return count_genes, count_exon, count_intron , float(len(genes_all)), chr

def process_feature(feature):
	feat=[]
	feat_all=[]
	chr=[]
	for line in open(infile,'r'):
		token = line.strip().split('\t')
		if token[2] == feature:
			if token[0] not in chr:
				feat.append([])
				chr.append(token[0])
				feat[len(feat)-1].append(abs(int(token[3])-int(token[4])))
				feat_all.append(abs(int(token[3])-int(token[4])))
			else:
				feat[len(feat)-1].append(abs(int(token[3])-int(token[4])))
				feat_all.append(abs(int(token[3])-int(token[4])))
				
	print "%s     %20s %30s %30s %30s %30s %30s" %(str(feature),str(len(feat_all)),str(np.sum(feat_all)),str(round(np.average(feat_all),2)),str(np.max(feat_all)),str(np.min(feat_all)),str(round(len(feat_all)/float(gene_counts),2)))
	count=[]
	for i in feat:
		count.append(len(i)) 
	return count

if __name__ == "__main__":
    
    ### get the options
    infile, gfile = options(sys.argv[1:])
    
    print infile
    
    ### read infile in memory
    inf = open(infile,"r")
    text = inf.read()
    inf.close()

    # search key words in the file
    # gene 
    print "Gene present in the GTF file:"+"\t"+str(search_file('gene'))
    # transcript
    print "Transcript present in the GTF file:"+"\t"+str(search_file('transcript'))
    # exon
    print "Exon present in the GTF file:      "+"\t"+str(search_file('exon'))
    # intron
    print "Intron present in the GTF file:    "+"\t"+str(search_file('intron'))
    #mRNA
    print "mRNA present in the GTF file:      "+"\t"+str(search_file('mRNA'))
	#UTR
    print "UTR present in the GTF file:      "+"\t"+str(search_file('UTR'))
    #FPKM
    print "FPKM present in the GTF file:      "+"\t"+str(search_file('FPKM'))
    #RPKM
    print "RPKM present in the GTF file:      "+"\t"+str(search_file('RPKM'))
    
    # process data
    genes,exon,intron, gene_counts, chr =process_data(infile)
    
    ### do the plotting
    col = []
    name = []
    import matplotlib.pyplot as plt
    N=len(exon)
    ind = np.arange(N)    # the x locations for the groups
    width = 0.3       # the width of the bars: can also be len(x) sequence
    p1 = plt.bar(ind+0.5,exon,width, color=(255/float(255),0/float(255),0/float(255)))
    col.append(p1[0])
    name.append('Exon')
    p2 = plt.bar(ind+0.5+width,intron,width, color=(255/float(255),241/float(255),183/float(255)))
    col.append(p2[0])
    name.append('Intron')
    def autolabel(rects):
    	# attach some text labels
    	for rect in rects:
    		height = rect.get_height()
    		plt.text(rect.get_x()+rect.get_width()/2., 0.5*height, str(round(100*height/float(sum(chr)),2))+'%',ha='center', va='bottom')

 
    
    # process features
    
    
    # gene
    if search_file('gene') == True:
    	gene_count=process_feature('gene')
    	p3 = plt.bar(ind+0.5+2*width,gene_count,width, color=(100/float(255),155/float(255),128/float(255)))
    	col.append(p3[0])
    	name.append('gene')
	
    # transcript
    if search_file('transcript') == True:
    	transcript=process_feature('transcript')
    	p3 = plt.bar(ind+0.5+2*width,transcript,width, color=(52/float(255),255/float(255),128/float(255)))
    	col.append(p3[0])
    	name.append('Transcript')
    # mRNA
    if search_file('mRNA') == True:
    	mRNA=process_feature('mRNA')
    	p4 = plt.bar(ind+0.5+3*width,mRNA,width, color=(233/float(255),138/float(255),255/float(255)))
    	col.append(p4[0])
    	name.append('mRNA')
    # UTR
    if search_file('UTR') == True:
    	UTR=process_feature('UTR')
    	p5 = plt.bar(ind+0.5+4*width,UTR,width, color=(62/float(255),248/float(255),255/float(255)))
    	col.append(p5[0])
    	name.append('UTR')
    text=''
    
        		
    #autolabel(p1)
    plt.xlabel('Chromosomes')
    plt.ylabel('Feature counts')
    plt.title('Feature count distrubution')
    plt.xticks(ind+width/2.+0.5, chr,rotation='10',fontsize=10 )
    plt.legend( col, name,bbox_to_anchor=(1.05, 1.05))
    plt.savefig(infile+".png")
    #plt.show()
    
### process the genome file