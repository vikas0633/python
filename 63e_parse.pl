#!/usr/bin/perl
use warnings;
use strict;
use diagnostics;

### declare variables
my $input1 = $ARGV[0];
my $input2 = $ARGV[1];
my $output1 = $ARGV[2];
my $output2 = $ARGV[3];
my $key; #key of entry in hash
my %entries; #hash
my @a1 = ();
my $row;

### open files
open (INFILE1, "<$input1") or die $!;
open (INFILE2, "<$input2") or die $!;
open (OUTFILE1, ">$output1") or die $!;
open (OUTFILE2, ">$output2") or die $!;


### populate a hash with data from INFILE1:
while(my $line1 = <INFILE1>) {
	chomp($line1);
	my @a = split "\t", $line1;
	my $sequence = $a[0];
	$entries{$sequence} = $line1;
}



#### extract from INFILE1 file based on INFILE2:
while (my $line2 = <INFILE2>) {
	chomp($line2);
	my @a = split "\t", $line2;
	$key = $a[10];
	# print $key;
	if (exists $entries{$key}) {
	$entries{$key} =~ s/\r//;
	$line2 =~ s/\r//;
	print  OUTFILE1 "$entries{$key}\t$line2\n" ;
	}
}

### close files
close INFILE1;
close INFILE2;
close OUTFILE1;

### make and sort the final IGV file
open (OUTFILE1, "<$output1") or die $!;
while (my $line3 = <OUTFILE1>) {
	chomp($line3);
	my @a = split "\t", $line3;
	push (@a1, [@a[33..35,0,15..26]]);
}

### sort records according to column 1 numeric
my @a1_sorted = sort{$a->[0] cmp $b->[0] ||
					$a->[1] <=> $b->[1]} @a1;

### print sorted records as a tab-separated file
print OUTFILE2 "chromosome\tstart\tstop\tsequence\twt_mock\twt_mlot\tnfr5_mock\tnfr5_mlot\tnfr1_mock\tnfr1_mlot\tsymrk_mock\tsymrk_mlot\tccamk_mock\tccamk_mlot\tcyclops_mock\tcyclops_mlot\n";
$" = "\t"; # sets the array separator to tabs!
for $row (@a1_sorted) {
print OUTFILE2 "@$row" , "\n";
}

#  <chromosome>	<start>	<stop>	<sequence>	<wt_mock>	<wt_mlot>	<nfr5_mock>	<nfr5_mlot>	<nfr1_mock>	<nfr1_mlot>	<symrk_mock>	<symrk_mlot>	<ccamk_mock>	<ccamk_mlot>	<cyclops_mock>	<cyclops_mlot>


## close files
close OUTFILE1;
close OUTFILE2;
