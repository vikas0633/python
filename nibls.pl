#!/usr/bin/perl
#####################################################################################################
#LocusPocus is a free script, it is provided with the hope that you will enjoy, you may freely redistribute it at will. We would be greatful if you would keep these acknowledgements with it. 
#
# Dan MacLean
# dan.maclean@sainsbury-laboratory.ac.uk
#
# This program is free academic software; academic and non-profit
# users may redistribute it freely.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
#
# This software is released under GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
# see included file GPL3.txt
#
#


###Dont forget you will need ... 
#####################################################################################################
# Boost::Graph 
#Copyright 2005 by David Burdick
# Available from http://search.cpan.org/~dburdick/Boost-Graph-1.2/Graph.pm
#Boost::Graph is free software; you can redistribute it and/or modify it under the same terms as Perl itself.
#####################################################################################################



use strict;
use warnings;
use Boost::Graph;
use Getopt::Long;


my $usage = "usage: $0 -f GFF_FILE [options]\n\n -m minimum inclusion distance (default 5)\n -c clustering coefficient (default 0.6) -b buffer between graphs (default 0)\n";

my $gff_file ;
my $min_inc = 5;
my $clus = 0.6;
my $buff = 0;

GetOptions(

	'c=f'          => \$clus,
	'm=i' => \$min_inc,
  'f|file=s'     => \$gff_file,
	'b=i' => \$buff
) ;



die $usage unless $gff_file;


my $starttime = time;
warn "started $starttime\n";

## load in data
my %molecules; # stores starts and ends of srnas
open GFF, "<$gff_file";

while (my $entry = <GFF>){

	chomp $entry;
	my @data = split(/\t/,$entry);
	
	#### modifying it read IGV format or Bed format 
	$molecules{$data[0]}{$data[1]}{$data[2]} = 1;

}

close GFF;

warn "Data loaded...\nBuilding graphs and finding loci\nPlease be patient, this can take a while...\n";




foreach my $chromosome (keys %molecules){
	my $g = new Boost::Graph(directed=>0);
	my @starts = keys(%{$molecules{$chromosome}} );
	@starts = sort {$a <=> $b} @starts;  

	while (my $srna_start = shift @starts){   ## work from left most sRNA to right most, add to graph if they close enough


		foreach my $srna_end (keys %{$molecules{$chromosome}{$srna_start}}){


			###use new graph if the next srna is too far away from this one.. 
			if(defined $starts[0] and $srna_end + $min_inc + $buff < $starts[0]){


				##dump the info from the old graph
				if (scalar(@{$g->get_nodes()}) > 2){

					my $cluster_coeff = get_cc($g);
					if ($cluster_coeff >= $clus){
					 dump_locus($g, $cluster_coeff);
					}
				}


				$g = new Boost::Graph(directed=>0);

			}

			foreach my $e (keys %{$molecules{$chromosome}{$srna_start}}){   ### extra bit because all loci with same start and different end overlap by definition. but are not collected by main search below

				unless ($e eq $srna_end){
					my $sn = $chromosome. ':' . $srna_start . ':' . $srna_end;   ## turn coordinate of sRNA inro a node name
					my $en = $chromosome. ':' . $srna_start . ':' . $e;
					$g->add_edge(node1=>"$sn", node2=>"$en", weight=>'1');
				} 

			}

			foreach my $start (@starts){  ##build graph of overlaps
				my $new = 0;
				last if $start - $min_inc > $srna_end;
				if ($start + $min_inc < $srna_end){

					my $start_node = $chromosome . ':' . $srna_start . ':' . $srna_end;
					foreach my $end (keys %{$molecules{$chromosome}{$start}}){

						my $end_node = $chromosome . ':' . $start . ':' . $end;
						$g->add_edge(node1=>"$start_node", node2=>"$end_node", weight=>'1');
					}

				}
			}
		}
	}
}

warn "Loci printed\nFinished\n";

my $endtime = time;

my $elapsed = $endtime - $starttime;

warn "Time elapsed = $elapsed s\n";	

#########################################################################################
sub get_cc{   ## do cluster coeff calculation. No useful method anyway so self implemented NB, this is an undirected graph so k is n(n-1)/2

	my $graph = shift;

	my @component = @{$graph->get_nodes()}; #number of nodes
	my @clustering_coefficients;

	foreach my $vertex (@component)
	{

		my @neighbours = @{$graph->neighbors($vertex)};

		my %edges_in_graph;

		my $n = @neighbours; #n = the number of neighbours
		my $k = ($n * ($n - 1))/2;	#k = total number of possible connections

		my $e= 0; #actual number of connections within sub-graph



		foreach my $neighbour (@neighbours)
		{
			foreach my $neighbour_2 (@neighbours)
			{
			my $edge1 = "$neighbour\t$neighbour_2";
			my $edge2 = "$neighbour_2\t$neighbour";
				unless (exists $edges_in_graph{$edge1} or exists $edges_in_graph{$edge2})
				{
					if ($graph->has_edge($neighbour, $neighbour_2) or $graph->has_edge($neighbour_2, $neighbour))
					{
						++$e;
						$edges_in_graph{$edge1}=1;
						$edges_in_graph{$edge2}=1;
					}
				}
			}
		}

		if ($k >= 1)
		{	
		my $c = $e / $k;
		push @clustering_coefficients, $c;
		}
		else {push @clustering_coefficients, '0';}
	}

	my $graph_n = scalar(@clustering_coefficients);
	my $graph_cc = 0;
	foreach my $cc (@clustering_coefficients){ 

		$graph_cc = $graph_cc + $cc;

	}
	$graph_cc = $graph_cc / $graph_n;

	return $graph_cc;
}

############################################################################################################

sub dump_locus{

my $g = shift;
my $cc = shift;
my $chr;
my $start = 1000000000000000000000000000000000000000000000; 
my $end = -1;

foreach my $node (@{$g->get_nodes()}){

	$node =~ m/^(\S+):(\d+):(\d+)$/;
	$chr = $1;
	$start = $2 if $2 < $start;
	$end = $3 if $3 > $end;


}
my $gff = $chr . "\t" . 'graph_locus_finder' . "\t" .'graph_locus' . "\t" . $start . "\t" . $end . "\t" . $cc . "\t" . '.' . "\t" . '.' . "\t" . '.';

print $gff, "\n";  

}