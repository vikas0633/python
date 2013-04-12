"""
my biopython functions
	
 v. 1.7  - 2011-Sep-28 - return an ordered list from readfasta
 v. 1.6  - 2010-May-18 - add space to lines of quality files
 v. 1.5  - 2009-Jul-01 - redirect diagnostic output to stderr
 v. 1.44 - 2009-Jun-10 - strip function removes spaces only from name
 v. 1.42 - 2009-May-21 - changed strip function to remove trailing spaces
 v. 1.41 - 2009-May-21 - added readfasta_nocheck for clean files
 v. 1.3  - 2009-May-16 - added removedupes option
 v. 1.2  - 2009-May-16 - added stringblock to return string instead of print
						fixed condition where it would fail if fasta had blank at top
 v. 1.1  - 2009-Mar-05.

readfasta
	given a fasta file name, return a dictionary of name:sequence
printblock
	print a string on lines of a certain width

stringblock
	build and return a string with lines of a certain width 
	(use for writing to file)

beroe <at } mac { dot> com

"""
# 
# removedupes applies to the sequences, not the names. those are removed 
# as long as appendseqs is false. otherwise they are appended
import sys

def readfasta(inputfilename, appendseqs = False, removedupes = True, inorder=False):
	dataset={}
	backseqs={}
	orderedlist=[]
	continueprocessing=True
	dupecount=0
	skipcount=0
	seqname=""
	skipnext=False;  firstseq=False;  nbrf=False; skipseq=False
	try: 
		inFile = file(inputfilename, 'rU') 
		alllines = inFile.readlines() # trying a different approach
	# this is easier for the \r detection, but for very large files,
	# it is probably better to use the original formulation
	# plus, does readlines() work in 2.4??

	except IOError:
		print >> sys.stderr ,"# Can't find the file %s. Are you in the right directory?" % inputfilename
		printtoscreen=False # false means it will save to the outfile
		printhtml=False
		savetofile=False
		alllines=[]
		continueprocessing = False
		return {}
		
	if continueprocessing:
		if alllines[0].find('\r')>0:
			print >> sys.stderr, ""
			print >> sys.stderr , "# MacFormat (use unix file when possible)."
			lines=alllines[0].split('\r')
			# print "Found %d lines" % len(lines)	

		else:
			lines=alllines

		# print len(lines)

		for inline in lines:
			line=inline.strip('\r').strip('\n')
			if skipnext: # needed for NBRF format
				skipnext=False
			else:
				if line and line[0]=="#":
					continue
				elif line and line[0]=='>':
					if firstseq:
						if dataset.get(seqname) in backseqs.keys(): # gets the sequence from the previous record n the loop
							dupecount+=1
							#print >> sys.stderr , "# Same sequence in",backseqs[dataset[seqname]],"and",seqname
							print >> sys.stderr, ".",
							if removedupes:
								del dataset[seqname]
						else:
							backseqs[dataset[seqname]]=seqname
					
					skipseq=False # if we are skipping a repeat, need to reset now
					if line.startswith('>DL;') or line.startswith('>P1;'):
						nbrf=True
						# print >> sys.stderr , 'NBRF format'
						seqname=line.split()[1].strip(' ')  # defaults to space
						skipnext=True
					else:
						seqname=line[1:].strip(' ')
						# print >> sys.stderr , seqname
					if seqname in dataset.keys() and (not appendseqs):
						print >> sys.stderr , "# Duplicate sequence name skipped" , seqname
						skipcount+=1
						skipseq=True # seqname already exists -- skip to next >	
					firstseq=True

				elif skipseq:
					continue

				elif firstseq: # we know it's not the first line
					# print >> sys.stderr ,  'we have sequence'
					try:
						dataset[seqname] +=  line
					except:
						dataset[seqname]  =  line
						orderedlist.append(seqname)
					# this should be modified to append
			# need to check the last sequence in the loop
		if dataset.get(seqname) in backseqs.keys(): # gets the sequence from the previous record n the loop
			dupecount+=1
			print >> sys.stderr , "# Same sequence in",backseqs[dataset[seqname]],"and",seqname
			if removedupes:
				del dataset[seqname]

		inFile.close()

		if nbrf: print >> sys.stderr , "# NBRF format" 
		if removedupes: print >> sys.stderr , "# Removed",dupecount,"identical sequences"
		if not appendseqs and skipcount>0: print >> sys.stderr , "# Removed",skipcount,"identical names"
		if inorder:
			return dataset,orderedlist
		else:
			return dataset



def printblock(seqstring, linelen = 60):
		
	"""print a chunk of text in lines of size linelen"""
	while seqstring:
		print seqstring[:linelen]
		seqstring=seqstring[linelen:]
	
def stringblock(seqstring, linelen = 60):
	outstring=''
	"""print a chunk of text in lines of size linelen"""
	while seqstring:
		outstring += seqstring[:linelen] +'\n'
		seqstring=seqstring[linelen:]
	return outstring

def readfasta_nocheck(inputfilename, appendseqs = False, removedupes = False, qualityfile=False, inorder=False):
	dataset={}
	orderedlist=[]
	continueprocessing=True
	seqname=""
	appender=""
	firstseq=False; 
	if qualityfile:
		appender=" "
	try: 
		inFile = file(inputfilename, 'rU') 
		alllines = inFile.readlines() # trying a different approach
	# this is easier for the \r detection, but for very large files,
	# it is probably better to use the original formulation
	# plus, does readlines() work in 2.4??

	except IOError:
		print >> sys.stderr , "# Can't find the file %s. Are you in the right directory?" % inputfilename
		printtoscreen=False # false means it will save to the outfile
		printhtml=False
		savetofile=False
		alllines=[]
		continueprocessing = False
		return {}
		
	if continueprocessing:
		if alllines[0].find('\r')>0:
			print >> sys.stderr ,""
			print >> sys.stderr , "# MacFormat (use unix file when possible)."
			lines=alllines[0].split('\r')
			# print "Found %d lines" % len(lines)	

		else:
			lines=alllines

		# print len(lines)

		for inline in lines:
			line=inline.strip('\r').strip('\n')
			if line and line[0]=='>':
				seqname=line[1:].strip(' ')
					# print >> sys.stderr , seqname
# should just use not(seqname=="")
			elif not(seqname==""): # we know it's not the first line
				# print >> sys.stderr ,  'we have sequence'
				try:
					dataset[seqname] +=  line + appender
				except:
					dataset[seqname]  =  line + appender
					orderedlist.append(seqname)
					

		inFile.close()

	if inorder:
		return dataset,orderedlist
	else:
		return dataset

