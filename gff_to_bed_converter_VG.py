import sys, re

assert sys.version_info[:2] >= ( 2, 4 )

def get_bed_line( chrom, name, strand, blocks ,last_t_st, last_t_en, last_ID, last_note, CDS):
	""" Returns a BED line for given data. """

	
	if len( blocks ) == 0:
		# Use simple BED format if there is only a single block:
		#   chrom, chromStart, chromEnd, name, score, strand
		#
		start, end = blocks[0]
		return "%s\t%i\t%i\t%s\t0\t%s\n" % ( chrom, start, end, name, strand )

	#
	# Build lists for transcript blocks' starts, sizes.
	#
	
	### get the block start and end
	CDS_min = sys.maxint
	CDS_max = -1
	for cds_st, cds_en in CDS:
		if cds_st < CDS_min:
			CDS_min = cds_st
		if cds_en > CDS_max:
			CDS_max = cds_en
	
	# Get transcript start, end.

	t_start = sys.maxint
	t_end = -1
	for block_start, block_end in blocks:
		if block_start < t_start:
			t_start = block_start
		if block_end > t_end:
			t_end = block_end

	# Get block starts, sizes.
	block_starts = []
	block_sizes = []
	for block_start, block_end in blocks:
		block_starts.append( str( block_start - t_start ) )
		block_sizes.append( str( block_end - block_start ) )
	
	
	#
	# Create BED entry.
	# Bed format: chrom, chromStart, chromEnd, name, score, strand, \
	#               thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
	#
	# Render complete feature with thick blocks. There's no clear way to do this unless
	# we analyze the block names, but making everything thick makes more sense than
	# making everything thin.
	#
	return "%s\t%i\t%i\t%s\t0\t%s\t%i\t%i\t0\t%i\t%s\t%s\t%s\t%s\n" % \
			( chrom, t_start, t_end, name, strand, CDS_min, CDS_max, len( block_starts ), ",".join( block_sizes ), ",".join( block_starts ), last_ID, last_note )

def __main__():
	input_name = sys.argv[1]
	output_name = sys.argv[2]
	skipped_lines = 0
	first_skipped_line = 0
	out = open( output_name, 'w' )
	i = 0
	cur_transcript_chrome = None
	cur_transcript_id = None
	cur_transcript_strand = None
	cur_transcripts_blocks = [] # (start, end) for each block.
	last_note = ' '
	CDS=[]
	last_t_st = -1
	last_t_en = -1
	last_ID = ''
	note = ' '
	for i, line in enumerate( file( input_name ) ):
		line = line.rstrip( '\r\n' )
		elems = line.split( '\t' )
		if line and not line.startswith( '#' ):
			if elems[2]!="mRNA":
				match = re.search(r'Parent=.+',line)
				if match:
					match = match.group().split(';')[0].replace('Parent=','')
					t_id = match
			
				# GFF format: chrom source, name, chromStart, chromEnd, score, strand, attributes
				start = str( long( elems[3] ) - 1 )
				coords = [ long( start ), long( elems[4] ) ]
				strand = elems[6]
				if strand not in ['+', '-']:
					strand = '+'
				''' original code
				#attributes = parse_gff_attributes( elems[8] )
				#t_id = attributes.get( "transcript_id", None )
				'''
				### following is modified to work with transcript_exon file from the RAP-DB
				### original code continues
				if not t_id:
					#
					# No transcript ID, so write last transcript and write current line as its own line.
					#
			
					# Write previous transcript.
					if cur_transcript_id:
						# Write BED entry.
						out.write( get_bed_line( cur_transcript_chrome, cur_transcript_id, cur_transcript_strand, cur_transcripts_blocks ) )
				
					# Replace any spaces in the name with underscores so UCSC will not complain.
					name = elems[2].replace(" ", "_")
					out.write( get_bed_line( elems[0], name, strand, [ coords ] ) )
					
				
				
				# There is a transcript ID, so process line at transcript level.
				if t_id == cur_transcript_id:
					# Line is element of transcript and will be a block in the BED entry.
					cur_transcripts_blocks.append( coords )
					last_t_st = t_st
					last_t_en = t_en
					last_ID = ID
					last_note = note
					if elems[2].strip()=="CDS": ###get the max min for thick labeling
						CDS.append( coords )
					continue
			
				
				#
				# Line is part of new transcript; write previous transcript and start
				# new transcript.
				#
				# Write previous transcript.
				
				if cur_transcript_id:
					# Write BED entry.
					out.write( get_bed_line( cur_transcript_chrome, cur_transcript_id, cur_transcript_strand, cur_transcripts_blocks ,last_t_st,\
					 last_t_en, last_ID, last_note, CDS) )
				
				# Start new transcript.
				cur_transcript_chrome = elems[0]
				cur_transcript_id = t_id
				cur_transcript_strand = strand
				cur_transcripts_blocks = []
				cur_transcripts_blocks.append( coords )  
				CDS=[]
				if elems[2].strip()=="CDS": ###get the max min for thick labeling
					CDS.append( coords )
										
			if elems[2]=="mRNA":
				t_st = ( long( elems[3] ) - 1 )
				t_en = ( long( elems[4] ) )
				match = re.search(r'ID=.+;',line)
				ID = match.group().split(';')[0].replace('ID=','')
				ID = ID.split('-')[0].replace('t','s')
				match = re.search(r'Note=.+;',line)
				if match:
					note = match.group().split(';')[0].replace('Note=','')
				if len(note) < 2:
					note = ' '
				continue
					
		else:
			skipped_lines += 1
			if not first_skipped_line:
				first_skipped_line = i + 1
    
	# Write last transcript.
	if cur_transcript_id:
		# Write BED entry.
		out.write( get_bed_line( cur_transcript_chrome, cur_transcript_id, cur_transcript_strand, cur_transcripts_blocks ,last_t_st, last_t_en, last_ID,\
		 last_note, CDS) )
	out.close()
	info_msg = "%i lines converted to BED.  " % ( i + 1 - skipped_lines )
	if skipped_lines > 0:
		info_msg += "Skipped %d blank/comment/invalid lines starting with line #%d." %( skipped_lines, first_skipped_line )
	print info_msg

if __name__ == "__main__": __main__()
