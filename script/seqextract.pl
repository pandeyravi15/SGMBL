#!/user/bin/perl
#print"enter the sequence file name:";
$fname=$ARGV[0];
open(FILE,$fname);

$fname1=$ARGV[1];
open(FILE1,$fname1);

$ofile = $ARGV[2];
open(out, ">$ofile");

#****************************************************
while( $line = <FILE> )
{
  chomp $line;
  $line =~ s/\n//;

 if( $line !~ /^>/) { $dna  .= $line; }

}

$s=length($dna);
print "$s\n";
#print out ">gene_cds_exon\n";
#print  out9 "$s\n";
#print out "$dna\n";


#*******************************************************************************
@D =<FILE1>;
$si=@D;

$a=0;$b=0;
$c=0;$d=0;
$e=0;$f=0;

foreach $D(@D)
{
	chomp $D;								
	@C = split(/\s+/,$D);

@arr1[$a++]=$C[0];
@arr2[$b++]=$C[1];
@arr3[$c++]=$C[2];
@arr4[$d++]=$C[3];
@arr5[$e++]=$C[4];
@arr6[$f++]=$C[5];
@arr7[$g++]=$C[6];
@arr8[$h++]=$C[7];
@arr9[$i++]=$C[8];
}

$len=$a;
for($i=0;$i< $len;$i++)
{
if($arr1[$i]!=$arr1[$i+1]){$QQ++;}


}
#****************************************************************************
print "$QQ\n";
#***************************************************************************
$k=0;
for($l=1;$l <= $QQ;$l++)
{
$k=0;print out ">$ARGV[0]_$l\n";
	for($i=0;$i< $len;$i++)
	{
		if($arr1[$i]==$l)
		{
			$f1=substr($dna,$arr2[$i],$arr5[$i]);
			print out "$f1\n";
			print "$arr1[$i]\t$arr2[$i]\t$k\t";
			$k+=length($f1);
			$j=$k-1;
			print "$j\n";
		}
	}

}
#print "$k\n";







