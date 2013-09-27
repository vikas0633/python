#-----------------------------------------------------------+
#                                                           |
# classCigar.py - process cigar string                      |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# AUTHOR: Vikas Gupta                                       |
# CONTACT: vikas0633@gmail.com                              |
# STARTED: 09/06/2013                                       |
# UPDATED: 09/06/2013                                       |
#                                                           |
# DESCRIPTION:                                              | 
# Short script to convert and copy the wheat BACs           |
# Run this in the parent dir that the HEX* dirs exist       |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+



def process(string):
    count = ''
    read_length = 0
    matches = 0
    insertions = 0
    deletions = 0
    
    for letter in string:
        if letter in 'MIDNSHPX=':
            
            ### count matches/mismatches
            if letter == 'M':
                matches += int(count) 
            ### count insertions
            if letter == 'I':
                insertions += int(count)
            ### count deletions
            if letter == 'D':
                deletions += int(count)   
            read_length += int(count)
            count = ''
            
        else:
            count += letter
    return read_length, matches, insertions, deletions

class Cigar:
    
    def __init__(self, string):
        self.string = string
        self.rl, self.m, self.i, self.d =  process(string)  
    
    def __str__(self):
        return self.string
    
    def read_length(self):
        return self.rl
    
    def matches(self):
        return self.m
    
    def insertions(self):
        return self.i
    
    def deletions(self):
        return self.d 