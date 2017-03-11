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
@arr4[$d++]=$C[3];
@arr5[$e++]=$C[4];
@arr6[$ee++]=$C[5];
@arr7[$dd++]=$C[6];	
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
@ar3[$q++]=$C1[2];
@ar4[$r++]=$C1[3];
@ar5[$z++]=$C1[4];
@ar6[$t++]=$C1[5];
@ar7[$u++]=$C1[6];
@ar8[$v++]=$C1[7];
}
$len2=$o;
#print "$len1\n";

#************************************************

for($i=0;$i<$len1;$i++)
{
$flag=0;
$j=0;
while($j<$len2 && $flag==0)
{

#if(($arr1[$i] == $ar1[$j])||($arr2[$i] == $ar2[$j])){                                            #compare both start and end position
if(($arr2[$i] eq $ar1[$j])){
#$ar3[$i]=$ar3[$i]+$arr2[$j];
#$ar4[$i]=$ar4[$i]+$arr2[$j];
#$ar5[$i]=$ar4[$i]-$ar3[$i];
#$arr6[$i]=$ar2[$j];
print out "$arr1[$i]\t$arr2[$i]\t$ar2[$j]\t$ar3[$j]\t$ar4[$j]\t$ar5[$j]\t$ar6[$j]\t$ar7[$j]\t$ar8[$j]\n";

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
#print "$w\t$w1\n";





