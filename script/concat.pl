#!/user/bin/perl
#print"enter the sequence file name:";
$fname=$ARGV[0];
open(FILE,$fname);

#$fname1=$ARGV[1];
#open(FILE1,$fname1);

$ofile = $ARGV[1];
open(out, ">$ofile");

$ofile1 = $ARGV[2];
open(out1, ">$ofile1");
#****************************************************
$s=0;$i=1;

while( $line = <FILE> )
{
  chomp $line;
  $line =~ s/\n//;
 #if( $line =~ /^>/) { $dna  = $line; print out "$dna\n";}
 if( $line !~ /^>/) { $dna  .= $line; }

   else
      {
	#print "$i\t$t\n";
	$s=length($dna);
	$len=$s-$t;
	if($i==1){$s=0;}
	print out "$i\t$s\n";
	$k[$i]=$s;
	 $i++;
      }
}
$len=$i;
#$DNA=reverse($dna);
$s=length($dna);
#print "$len\n";
#print out ">gene_cds_exon\n";
#print  out9 "$s\n";
#print out "$dna\n";

for($i=1;$i<($len-1);$i++)
{print out1"$k[$i]\t$k[$i+1]\n";}
print out1"$k[$len-1]\t$s\n";
