$fname = $ARGV[0];
open (FILE, $fname);

$top = $ARGV[1];
$ofile = $ARGV[2];
open(out, ">$ofile");



while( $line = <FILE> )
{
  chomp $line;
   @fields = split (/\t/, $line); 
#print "@fields\n";

@arr1[$i++] = $fields[0]; 
@arr2[$j++] = $fields[1]; 
@arr3[$k++] = $fields[2];
@arr4[$l++] = $fields[3]; 
@arr5[$m++] = $fields[4]; 
@arr6[$n++] = $fields[5];

}
$length = $i; 
#print "$length\n";
$b[0]=0;
for($i = 0; $i < $length; $i++)
{
	if($arr1[$i]!=$arr1[$i+1])
	{
		$w++;$donar++;
		$a[$donar]=$b[$donar-1];
		$b[$donar] =$w;
		print "$donar\t$a[$donar]\t$b[$donar]\n";
		
	}
	else
	{
		#print  out "$a[$i]\t$b[$i]\t$c[$i]\n";
		$w++;
	}

}

print "$donar\n";
###################################################################################################################################

$j=0;
for($r = 1; $r <= $donar; $r++)
{
	$k=0;	

  	for($i = $a[$r]; $i < $b[$r]; $i++)
  	{

		#if(($arr1[$i] == $r))	 		
	 	 # {	
			$u[$k]=$arr1[$i];$v[$k]=$arr2[$i];$w[$k]=$arr3[$i];
			$k++;
		 		
			
	 	  # }

                      
	 			
		}
				
	for($k = 0; $k < $top; $k++)
 	{
		print  out "$u[$k]\t$v[$k]\t$w[$k]\n";
	}			
		
}
		
	
#print  "$count\n";	
###################################################################################################################################	
#sort -k3,3 kk2 >kk3
 #sort -k1,1 -k3,3 kk3 >kk4
###################################################################################################################################			











