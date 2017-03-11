my $f = $ARGV[0];
open(FILE,$f);

$ofile = $ARGV[2];
open(out, ">$ofile");



@D =<FILE>;
$si=@D;
#print "$D[2]";
$COUNT = 0;
$i=0;

foreach $D(@D)
{
	chomp $D;						
	@C = split(/\s+/,$D);
@arr1[$a++]=$C[0];
@arr2[$b++]=$C[1];
@arr3[$c++]=$C[2];

	
}
$len1=$a;
#**************************************************
my $f1 = $ARGV[1];
open(FILE1,$f1);

@D1 =<FILE1>;
$s=@D1;
#print "$D[2]";
$COUNT1 = 0;

foreach $D1(@D1)
{
	chomp $D1;						
	@C1 = split(/\s+/,$D1);
@ar1[$o++]=$C1[0];
@ar2[$p++]=$C1[1];
@ar3[$pp++]=$C1[2];
}
$len2=$o;
print "$len1\n";

#************************************************

for($i=0;$i<$len1;$i++)
{
$flag=0;
$j=0;
while($j<$len2 && $flag==0)
{

#if(($arr1[$i] == $ar1[$j])||($arr2[$i] == $ar2[$j])){                                            #compare both start and end position
if(($arr1[$i] == $ar1[$j])){
#$ar3[$i]=$ar3[$i]+$arr2[$j];
#$ar4[$i]=$ar4[$i]+$arr2[$j];
#$ar5[$i]=$ar4[$i]-$ar3[$i];
#$arr6[$i]=$ar2[$j];									  #compare only end position
print out "$arr1[$i]\t$ar2[$j]\t1\t$arr2[$i]\t$arr3[$i]\n";
$w++;
#print "$ar1[$i]\t$ar2[$i]\t$ar3[$i]\t$ar4[$i]\t$ar5[$i]\t$ar6[$i]\t$ar7[$i]\t$ar8[$i]\t$ar9[$i]\t$arr2[$j]\n";
$flag=1
}

$j++;
}
if($flag==0)
{
#print out "$arr1[$i]\t$arr2[$i]\tFALSE\n";
$w1++;
}
}#
print "$w\t$w1\n";





