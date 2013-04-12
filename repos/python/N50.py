#!/usr/bin/python
import sys

def N50(numlist):
    """
    Abstract: Returns the N50 value of the passed list of numbers. 
    Usage:    N50(numlist)

    Based on the Broad Institute definition:
    https://www.broad.harvard.edu/crd/wiki/index.php/N50
    """
    numlist.sort()
    newlist = []
    for x in numlist :
        newlist += [x]*x
    # take the mean of the two middle elements if there are an even number
    # of elements.  otherwise, take the middle element
    if len(newlist) % 2 == 0:
        medianpos = len(newlist)/2  
        return float(newlist[medianpos] + newlist[medianpos-1]) /2
    else:
        medianpos = len(newlist)/2
        return newlist[medianpos]

assert N50([2, 2, 2, 3, 3, 4, 8, 8]) == 6

lengths = []
for line in sys.stdin :
    lengths.append(int(line))
print N50(lengths)