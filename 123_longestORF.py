
import sys, re
def orf_find(st0):

    seq_0=""
    for i in range(0,len(st0),3):
        if len(st0[i:i+3])==3:
            seq_0 = seq_0 + st0[i:i+3] + ' '

    ms_1 =[m.start() for m in re.finditer('ATG', seq_0)]
    ms_2 =[m.start() for m in re.finditer('(TAA)|(TAG)|(TGA)', seq_0)]

    def get_next(arr,value):
        for a in arr:
            if a > value:
                return a
        return -1




    codons = []
    start_codon=ms_1[0]
    while (True):
        stop_codon = get_next(ms_2,start_codon)
        if stop_codon == -1:
            break
        codons.append((start_codon,stop_codon))
        start_codon = get_next(ms_1,stop_codon)
        if start_codon==-1:
            break

    max_val = 0
    selected_tupple = ()
    for i in codons:
        k=i[1]-i[0]
        if k > max_val:
            max_val = k
            selected_tupple = i

    print "selected tupple is ", selected_tupple

    final_seq=seq_0[selected_tupple[0]:selected_tupple[1]+3]
    print "The longest orf length is " + str(max_val)
    
    return final_seq.replace(' ','')
    



output_file = open('Longorf.txt','w')
for line in open(sys.argv[1]):
    if len(line)>0 and not line.startswith('#'):
        if line.startswith('>'):
            output_file.write(str(line))
        else:
            output_file.write(str(orf_find(line)))


output_file.close()
