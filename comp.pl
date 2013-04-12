open(INFILE,"/Users/vikas0633/Desktop/plant/2011_week49/jorge/srna_data.txt") or die "Coudln't access target file\n";
open(OUTFILE,">/Users/vikas0633/Desktop/plant/2011_week49/jorge/summaryup") or die "Coudln't access target file\n";
@nuc=("A","T","C","G");
while(<INFILE>)
{
	if(/\t\\ /)
	{
	chomp;
	@array=split(/\t/);
	print("$array[2]\n");
	$pal=substr $array[0],0,3,"";		
	foreach(@nuc)
		{
			$cont{$_}+=($array[0]=~s/$_//g);
			$cont2{$_}+=($pal=~s/$_//g);
			print("$cont{$_}\t$cont2{$_}")
		}
	}
}
foreach(@nuc)
{
	$tam+=$cont{$_};
	$tam2+=$cont2{$_};
}
foreach(@nuc)
{
#$num2=$cont{$_}/$tam;
#$num=$cont2{$_}/$tam2;
print("$_\t$num\t$num2\n");
}
close(OUTFILE);