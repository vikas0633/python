## Update: 20110204
## Auther: vikas gupta
## using variable window = no of tasi-RNAs


#!/usr/bin/perl
# This perl script is developed to identify small RNAs clustered in 21-nt increments.
# A 231-bp fragment with 2-nt 3'-overhang downstream from the 5' start site of each small RNA signature was used for the
# calculation of the number of small RNAs present (n),and the number of small RNAs in 21-nt increments relative to the
# start site (k). The input file (sRNAmapping.txt) should contain the start position of each small RNA signature in
# terms of chromosome, coordinate, and strand (1: Watson; -1: Crick) which are separated by Tab as shown below.
# -----------------------------------------------
# chr1	11927	-1	TGGCGATGATGATCAAT
# chr1	28152	1	CATCCTTCGATGTTGTG
# chr1	42408	1	CTCTTAGCTAAGAGCCA
#------------------------------------------------
# Then the p-value (p) of obtaining equal or more than k phased small RNAs is calculated based on hypergeometric
# distribution (see methods).
# Sep. 2006

open(data, "sRNAmapping.txt") || die "Cannot open sRNAmapping.txt";
open(out, ">all_sRNA21nt_out.txt"); # outfile of all small RNA clusters
open(out2, ">p0.001_sRNA21nt_out.txt"); # outfile of small RNA clusters with p<0.001

while($line=<data>) # input small RNA mapping data
{
    chomp $line;
    ($chr, $cor, $str, $seq)=split(/\t/, $line);
    $sp=$chr.",".$cor.",".$str;
    $sh{$sp}=$line;
}

foreach $ccs(keys%sh)
{	
	for($window =5; $window<20;$window++)
	{
    	$n=-1;
    	$k=-1;
    	($chr, $cor, $str)=split(/,/, $ccs);
    	if ($str==1) # small RNAs on the Watson strand
    	{
        	$ss=$cor;
        	$ee=$cor+ ($window)*21 -1;
        	for ($i=$ss; $i<=$ee; $i++) # calculate n on the sense strand
        	{
            $x=$chr.",".$i.","."1";
            if ($sh{$x}=~ m/\S+/)
                {
                    $n=$n+1;
                }
     	   }
        	for ($i=$ss; $i<=$ee; $i=$i+21) # calculate k on the sense strand
        	{
          	  	$x=$chr.",".$i.","."1";
            	if ($sh{$x}=~ m/\S+/)
            	{
                	$k=$k+1;
            	}
        	}
        	for ($j=$ss-2; $j<=$ee-2; $j++) # calculate n on the antisense strand
        	{
            	$x=$chr.",".$j.","."-1";
            	if ($sh{$x}=~ m/\S+/)
            	{
                	$n=$n+1;
            	}
        	}
        	for ($j=$ss+18; $j<=$ee-2; $j=$j+21) # calculate k on the antisense strand
        	{
            	$x=$chr.",".$j.","."-1";
            	if ($sh{$x}=~ m/\S+/)
            	{
                	$k=$k+1;
            	}
        	}
    	}
    	elsif ($str==-1) # small RNAs on the Crick strand   
    	{
        	$ee=$cor;
        	$ss=$cor-($window*21)+1;
        	for ($i=$ss; $i<=$ee; $i++) # calculate n on the sense strand
        	{
            	$x=$chr.",".$i.","."-1";
            	if ($sh{$x}=~ m/\S+/)
            	{
                	$n=$n+1;
            	}
        	}
        	for ($i=$ss+20; $i<=$ee; $i=$i+21) # calculate k on the sense strand
        	{
            	$x=$chr.",".$i.","."-1";
            	if ($sh{$x}=~ m/\S+/)
            	{
                	$k=$k+1;
                }
        	}
        	for ($j=$ss+2; $j<=$ee+2; $j++) # calculate n on the antisense strand
        	{
            	$x=$chr.",".$j.","."1";
            	if ($sh{$x}=~ m/\S+/)
            	{
                	$n=$n+1;
            	}
        	}
        	for ($j=$ss+2; $j<=$ee+2; $j=$j+21) # calculate k on the antisense strand
        	{
            	$x=$chr.",".$j.","."1";
            	if ($sh{$x}=~ m/\S+/)
            	{
                	$k=$k+1;
            	}
        	}
    	}
    	$p=0;
    	for ($w=$k; $w<=(2*$window-1); $w++) # calculate p-value from n and k
    	{
        	$c=1;
        	$rr=1;
        	$rw=1;
        	for ($j=0; $j<=$w-1; $j++)
        	{
            	$c=$c*($n-$j)/($j+1);
        	}
        	for ($x=0; $x<=$w-1; $x++)
        	{
            	$rr=$rr*(21-$x)/((2*21*$window)-1-$x);
        	}
        	for ($y=0; $y<=$n-$w-1; $y++)
        	{
            	$rw=$rw*((2*21*$window)-22-$y)/((2*21*$window)-1-$w-$y);
        	}
        	$pr=$c*$rr*$rw;
        	$p=$p+$pr;
    	}

    	print out "$sh{$ccs}\t$n\t$k\t$window\t$p\n"; #output all small RNA clusters with their p values

    	if ($p<0.001) #select and output small RNA clusters with p<0.001
    	{
        	print out2 "$sh{$ccs}\t$n\t$k\t$window\t$p\n";
    	}
	}
}
