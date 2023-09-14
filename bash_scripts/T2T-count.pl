$file=shift;
$chr_file=shift;
$length=shift;
$cutoff=shift;

die "Usage: $0 <tidk_result> <fa.fai> <length> >output
	<tidk_result> -- telemere results obtained from tidk find
	<fa.fai> -- fasta index file for the contig/scaffold
	<length> -- length to check at begining and end. [default: 100KB]
	<repeat_count> -- Number of repeats within <length> to considered as telomere contig [default: 100]

	With default values, any contig which contains at least 100 repeat units at 100kb start and end of sequence will be considered as T2T\n\n" if (!$file || !$chr_file);

$length=100000 if (!$length);
$cutoff=100 if (!$cutoff);

print STDERR "Analysis is proceeding with\n\t\ttelomer file-- $file\n\t\tchromosome length file-- $chr_file\n\t\tCutoff lenght-- $length\n\n";

($name, $haplot)=$chr_file=~/(\w+?)\_(hap\d+)/;

open KK, $chr_file or die "$chr_file --$!\n";

open MM, ">>Plot_T2T.txt" or die "Cannot write files in the given location - $!\n";

foreach (<KK>)
{
  chomp;
  @temp=split(/\t/, $_);
  $hash{$temp[0]}=$temp[1]-$length;
  $total_seq++;
}

close KK;

open KK, $file or die "$file --$!\n";
<KK>;
while (<KK>)
{
  chomp;
  @temp=split(/,/, $_);
  $input{$temp[0]}=1;
  if ($temp[1]<=$length)
  {
	$start{$temp[0]} += $temp[2]+$temp[3];
  }
  if ($temp[1]>=$hash{$temp[0]})
  {
	$end{$temp[0]}+= $temp[2]+$temp[3];
  }
}
close KK;

$t2t=$p5=$p3=$notelo=0;
#$t2tc=$p5c=$p3c=0;

print "#Name\tLength\t5'_repeatCount\t3'_repeatCount\n";

foreach (keys %input)
{
	#print $_,"\t",$start{$_},"\t",$end{$_},"\n";
	if ($start{$_} >= $cutoff && $end{$_} >= $cutoff)
	{
		print $_,"\t",$hash{$_}+$length,"\t", $start{$_},"\t",$end{$_},"\n";
		$t2t++;
	}
	if ($start{$_} >= $cutoff && $end{$_} < $cutoff)
	{
		print $_,"\t",$hash{$_}+$length,"\t", $start{$_},"\t",$end{$_},"\n";
		$p5++;
	}
	if ($start{$_} < $cutoff && $end{$_} >= $cutoff)
	{
		print $_,"\t",$hash{$_}+$length,"\t", $start{$_},"\t",$end{$_},"\n";
		$p3++;
	}
	if ($start{$_} < $cutoff && $end{$_} < $cutoff)
	{
		$notelo++;
	}
}

#print STDERR "Total sequences\t$total_seq\nT2T sequences\t$t2t\nTelomere at 5'\t$p5\nTelomere at 3'\t$p3\nNo telomere motif\t$notelo\n";

print MM "$name\t$haplot\tTotal\t$total_seq\n$name\t$haplot\tT2T\t$t2t\n$name\t$haplot\t5p\t$p5\n$name\t$haplot\t3p\t$p3\n$name\t$haplot\tNo_telo\t$notelo\n";
