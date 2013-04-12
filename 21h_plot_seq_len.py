### 21h_plot_seq_len.py - take a fasta file and plot sequence length -  /Users/vikas0633/Desktop/script/python


'''

'''
import sys

from Bio import SeqIO
sizes = [len(rec) for rec in SeqIO.parse(sys.argv[1], "fasta")]

print "Total number of sequences in the file:", len(sizes)
print "Total residue in the file:", sum(sizes)

import pylab
pylab.hist(sizes, bins=100, color='green')
pylab.title("%i fasta reads\nLengths %i to %i" \
            % (len(sizes),min(sizes),max(sizes)))
pylab.xlabel("Read length (bp)")
pylab.ylabel("Count")

import datetime
now = datetime.datetime.now()
pylab.savefig(str(now.strftime("%Y-%m-%d_%H%M_"))+'.png')



'''
>>> from Bio import SeqIO
>>> help(SeqIO.parse)

    
     - handle   - handle to the file, or the filename as a string
                  (note older verions of Biopython only took a handle).
     - format   - lower case string describing the file format.
     - alphabet - optional Alphabet object, useful when the sequence type
                  cannot be automatically inferred from the file itself
                  (e.g. format="fasta" or "tab")
    
    Typical usage, opening a file to read in, and looping over the record(s):
    
    >>> from Bio import SeqIO
    >>> filename = "Fasta/sweetpea.nu"
    >>> for record in SeqIO.parse(filename, "fasta"):
    ...    print "ID", record.id
    ...    print "Sequence length", len(record)
    ...    print "Sequence alphabet", record.seq.alphabet
    ID gi|3176602|gb|U78617.1|LOU78617
    Sequence length 309
    Sequence alphabet SingleLetterAlphabet()
    
    For file formats like FASTA where the alphabet cannot be determined, it
    may be useful to specify the alphabet explicitly:
    
    >>> from Bio import SeqIO
    >>> from Bio.Alphabet import generic_dna
    >>> filename = "Fasta/sweetpea.nu"
    >>> for record in SeqIO.parse(filename, "fasta", generic_dna):
    ...    print "ID", record.id
    ...    print "Sequence length", len(record)
    ...    print "Sequence alphabet", record.seq.alphabet
    ID gi|3176602|gb|U78617.1|LOU78617
    Sequence length 309
    Sequence alphabet DNAAlphabet()
    
    If you have a string 'data' containing the file contents, you must
    first turn this into a handle in order to parse it:
    
    >>> data = ">Alpha\nACCGGATGTA\n>Beta\nAGGCTCGGTTA\n"
    >>> from Bio import SeqIO
    >>> from StringIO import StringIO
    >>> for record in SeqIO.parse(StringIO(data), "fasta"):
    ...     print record.id, record.seq
    Alpha ACCGGATGTA
    Beta AGGCTCGGTTA
    
    Use the Bio.SeqIO.read(...) function when you expect a single record
    only.
'''