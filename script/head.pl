#!/user/bin/perl
#print"enter the sequence file name:";
$fname=$ARGV[0];
open(FILE,$fname);

#$fname1=$ARGV[1];
#open(FILE1,$fname1);

$ofile = $ARGV[1];
open(out, ">$ofile");

#****************************************************
$s=0;
$i=1;
$w=1;
while( $line = <FILE> )
{
  chomp $line;
  $line =~ s/\n//;
 if( $line =~ /^>/) { $dna  = $line; print out "$w\t$dna\n";$w++;}
 #if( $line !~ /^>/) { $dna  .= $line; }

   else
      {
	#print "$i\t$t\n";
	$s=length($dna);
	if($s>0){$s==$s;}else{$s=0;}
	#print out"$i\t$s\n";
	$i++;
      }
}

#$DNA=reverse($dna);
$si=length($dna);
print "$si\n";
#print out ">gene_cds_exon\n";
#print  out9 "$s\n";
#print out "$dna\n";

