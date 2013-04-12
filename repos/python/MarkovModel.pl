@nuc=("A","C","G","T");
$limit=10;
$min_success=3;
$alpha=.001;
sub counter
{
#This sub counts has as inputs: array_reference, kmer_size. This will count the counts of each kmer of size kmer_size and return it as a hash, where the key is the kmer and its contained
#value is the count number.
	my @list=@{$_[0]};
	my $kmer=$_[1];
	my $tam=0;
	my %cont;
	foreach(@list)
	{
		$tam=length($_)-$kmer+1;
		if($tam<1)
		{
			print("Error, sequence smaller than $kmer-mer\n");
			return(0);
		} 
		for($i=0;$i<$tam;$i++)
		{
			$pal=substr $_,$i,$kmer;
			$cont{$pal}++;		
		}
	
	}
	return(%cont);
}
sub totals
{
#This sub counts has as inputs: hash_reference. This will count the appeareances of (k-1)-mers, where k is the length of the key and return it as a hash, where the key is the kmer and its contained
#value is the count number.For example, if you have a hash with all the 3-mers, it will return you a hash with all the counts of 2 mers.
	my %list=%{$_[0]};
	my $sum=0;
	if($list{A})
	{
		for my $key (keys %list)
		{
			$sum+=$list{$key};
		}
		$list{1}=$sum;
		return(%list);
	}
	my $sum=0;
	my $tmp;
	my %tam_lis;
	for my $key ( keys %list ) 
	{
        $base=substr($key,0,-1);
       	$tam_lis{$base}+=$list{$key};       	
	}
    return(%tam_lis)
}
sub markov
{
#This sub counts has as inputs: hash_reference1,hash_reference2. This will calculate the instance of each k-mer in hash_reference divided by the (k-1)-mers from hash reference 2, and will return a hash reference where the key is the kmer and its contained
#value is the count number. For example, will divide the instances of GAA by the number of GA. It gives the markov model of order-k frequences.
	my %counts=%{$_[0]};
	my %tam=%{$_[1]};
	if ($tam{1})
	{
		print("size: $tam{1}\n");
		for $key ( keys %counts )
		{
			$counts{$key}=$counts{$key}/$tam{1};
		}	
		return(%counts)
	}
	my $tmp=0;
	my $key="";
	$number = keys %counts;
	$number2=keys %tam;
	for $key ( keys %counts )
	{
		$base=substr($key,0,-1);
		$mm{$key}=$counts{$key}/$tam{$base};
	}
	return(%mm);
		
	
}
sub columns
{
	#This sub takes as parameters array_ref(each element of the array contains a sequence), NUM . It will calculate the NUM-order markov probability for every NUM-mer of every column, starting from the column NUM.
	#for example, if num is 3 and your sequences have length 5, it will calculate the probability of seing each triplete given the abundance of every 2 duplete, for every column starting from 1, until 3
	my @seq=@{$_[0]};
	my $word_size=$_[1];
	my $k,$n,$words,$data,$len,%col,%col2,%col3;	
	$len=length($seq[0])-$word_size+1;
	print("$seq[0]\t$len\n");
	if($len>0)
	{
		for($k=0;$k<$len;$k++)
		{
			print("Calculatin column $k..\n");
			@temp_a=();
			$testo=@seq;
			print("init:$testo\n");
			for($n=0;$n<@seq;$n++)
			{
				$temp_a[$n]=substr $seq[$n],$k,$word_size;
			}
			$tami=@temp_a;
			print("t:$tami\n");
			%col=counter(\@temp_a,$word_size);
			%col2=totals(\%col);
			for my $key (keys %col2)
			{
				print("$key\t$col2{$key}\n");	
			}
			%col3=markov(\%col,\%col2);
			$words='';
			for $key (keys %col3)
			{
				$words.="$key\t";
				$data.="$col3{$key}\t";
			}
			$data.="\n";
		}
		print("$words\n$data");
		return(0);
	}
	else
	{
		print"FAIL! kmer bigger than seq\n";
		return(0);
	}
	
}
sub mm_file
{
#This is me being lazy, takes an array_reference,NUM and caclulates de NUM-order Markov Model frequencies for the strings contained in that array.
	my @seq=@{$_[0]};
	my $word_size=$_[1];
	my %col=counter(\@seq,$word_size);
	my %col2=totals(\%col);
	for my $key (sort keys %col2)
	{
		print("$key\t$col2{$key}\n");	
	}
	print("\n");
	for my $key (sort keys %col)
	{
		print("$key\t$col{$key}\n");	
	}
	my %col3=markov(\%col,\%col2);
	print("\n");
	for my $key (sort keys %col3)
	{
		print("$key\t$col3{$key}\n");	
	}
	return(%col);
		
}
sub kmer
{
	#this sub takes as a parameter array_reference,NUM. It separates every string contained in the array in NUM-mer and (K-NUM)-mer, where K is the length of the sequence.
	my @seq=@{$_[0]};
	my $word_size=$_[1];
	my @mer,my @rest;
	my $i=0;
	foreach(@seq)
	{
		$mer[$i]=substr $_,0,$word_size,"";
		$rest[$i]=$_;
		$i++;
	}
	return(\@mer,\@rest);
}
sub complicated
{
	#This sub takes an array_ref full of sequences, and a kmer size, and a sequence size. It will count the nucleotide frequency in the sequence of sequence size, starting from the end.
	#It will return a hash that has k-mer's as keys, and nucleotides as second key. The reson for doing this is that you can build a contigency table with this, and then do a Chi test
	#or a G test, and see if the last nucleotides in your sequences are biased by the starting k-mer.For example, if the fact that your sequences Start with GAA, enriches for GAA en in the
	# whole sequences (proving if a sequence that starts with GAA have more GAA than expected)
	my @seq=@{$_[0]};
	my $word_size=$_[1];
	my $seq_size=$_[2];
	my $kmer=$_[3];
	my $key,%mer,@rest,%tmp,$key2;
	$seq_size=0-$seq_size;
	if(length($seq[0])<$word_size+$seq_size+1)
	{
		print("Error, the word size+the seq size are bigger than the actual seq\n");
		return(0);
	}
	foreach(@seq)
	{
		$key=substr $_,0,$word_size;
		@rest=();
		$rest[0]=substr $_,$seq_size;
		print("$rest[0]\n");	
		%tmp=();
		%tmp=counter(\@rest,$kmer);
		for $key2 (sort keys %tmp)
		{
			$mer{$key}{$key2}+=$tmp{$key2};
		}
		print("\n");
	}
	for $key (sort keys %mer)
	{
		for $key2 (sort keys %{$mer{$key}})
		{
			print("$mer{$key}{$key2}\t");
		}
		print("\n")
	}
	return(1);
}
sub binomial
{
	my $success=$_[0];
	my $total=$_[1];
	my $prob=$_[2];
	open(OUTFILE,">myScript.R");
	print OUTFILE ("num<-pbinom($success,$total,$prob,lower.tail=FALSE)\nwrite(num, file='tmp_R', append = FALSE)");	
	close(OUTFILE);
	`R --vanilla <myScript.R`;
	open(INFILE,"tmp_R") or die "Coudln't access target file\n";
	$_=<INFILE>;
	chomp;
	close(INFILE);
	return($_)
}
sub seq_sort
{
	#this sub take as an input an array reference with sequences in it, and a string, and will return all sequences that start with that sequence. For example
	#if the second parameter is GAA, it will return all sequences that start with GAA in the form of an array.
	my @array=@{$_[0]};
	my @array2;
	my $key=$_[1];
	foreach(@array)
	{
		if(/^$key/)
		{
			push(@array2,$_);
		}
	}
	return(@array2);
}
sub kmer_discovery
{
	my @seq=@{$_[0]};
	my $tam=@seq;
	if($tam<$limit)
	{
		print("Not enough sample size\n");
		return(0);
	}
	my $level=$_[1];
	my $ref1,my $ref2,my $total,my $tmp,my @seq1,my $prob;
	if($level)
	{
		($ref1,$ref2)=kmer(\@seq,1);
		my @seq1=@{$ref1};
		my %hash=counter(\@seq1,1);
		for my $key (keys %hash)
		{
			$total+=$hash{$key};
		}
		for $key (keys %hash)
		{
			if($hash{$key}>$min_success)
			{
				$tmp=$hash{$key}-1;
				$prob=binomial($tmp,$total,$prob{$key});
				if($prob<$alpha)
				{
					my @array2=seq_sort(\@seq,$key);
					($ref1,$ref2)=kmer(\@array2,1);
					@array2=@{$ref2};
					$tam=$level+1;
					for($i=0;$i<$level;$i++)
				{
					print("\t");
				}
				print("$key\t$prob\n");
					kmer_discovery(\@array2,$tam);
				}				
			}
			else
			{
				print("Not enough successes:$hash{$key}\n");
				return(0);
			}	
		}
	}
	else
	{
		foreach(@nuc)
		{
			print("$_\n");
			my @array2=seq_sort(\@seq,$_);
			($ref1,$ref2)=kmer(\@array2,1);
			@array2=@{$ref2};
			$tam=$level+1;
			kmer_discovery(\@array2,$tam);
		}
	}
}
sub get_prob
{
	#This reads the file prob from mt_srna folder, and puts in a hash where teh keys are the first columns,
	my %hash;
	open(INFILE,"prob") or die "Coudln't access target file\n";
	while(<INFILE>)
	{
		chomp;
		my @array=split /\t/;
		$hash{$array[0]}=$array[1];
	}
	close(INFILE);
	return(%hash)
}
$file=$ARGV[0];
print("$file\n");
#To do the kmer discovery you need to have some probabilities a priori.
%prob=get_prob;
open(INFILE,$file) or die "Coudln't access target file\n";
$i=0;
while(<INFILE>)
{
	if(/[ATCG]/)
	{
	chomp;
	$array[$i]=$_;
	$i++;
	}
}
print("Size\n$i\n");
close(INFILE);
#This calculates the single nucleotide frequencies in the first 3 columns.
($ref1,$ref2)=kmer(\@array,3);
@a_ref1=@{$ref1};
@a_ref2=@{$ref2};
mm_file(\@a_ref1,3);
#$numero=@a_ref2;
@array2=@{$ref2};
@array1=@{$ref1};
#%hash=counter($ref1,3);
#%hash2=counter($ref2,3);
#for $key (sort keys %hash)
#{
#	print("$key\t$hash{$key}\t$hash2{$key}\n")
#}
#This calculates the second order markov model frequencies for each sequence group of sequences starting with each nucleotide
#complicated(\@array,1,10,2);
#kmer_discovery(\@array,0)