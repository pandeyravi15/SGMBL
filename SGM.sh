#!/bin/bash
FILE=$1

F=$2

python script/revcomp.py $FILE >$FILE.rev  #making reverse complement of input file
if [ -f $F ] ; then
    rm $F
fi
start_time=$(date +%s)
perl script/concat.pl $FILE $FILE.pos.txt $FILE.loc.txt
echo "starting scoring read sets using SGM["$(date)"]"

perl script/head.pl $FILE $FILE.head2.txt
sed -i 's/>//g' $FILE.head2.txt
awk  '{gsub(".fna_","\t",$0); print;}' $FILE.head2.txt >$FILE.head.txt

for i in 1 2 3 4 5 6 7 8
do
	cc script/model5.c -lm
	for file in sgm_models/*.$i.param.txt; do 
	if [ -f "$file" ];then
   		./a.out "$file" $FILE $i $F.fw $FILE.loc.txt >out;
		./a.out "$file" $FILE.rev $i $F.rev $FILE.loc.txt >out;
	else
   		echo "File $file does not exist" 
	fi
done 
done
#cc model5.c -lm
#for file in models/*.1.param.txt; do ./a.out "$file" $FILE 1 $F $FILE.loc.txt;done >out


finish_time=$(date +%s)
echo "Done scoring read sets using SGM["$(date)"]"
echo "Time duration of scoring reads: $((finish_time - start_time)) secs."
sed -i 's/\///g' $F.fw
sed -i 's/sgm_models//g' $F.fw
sed -i 's/\///g' $F.rev
sed -i 's/sgm_models//g' $F.rev
#sed -i 's/model_fixedcluster//g' $F
cat $F.fw $F.rev >$F
sort -k1,1n $F >$F.2.txt




sort -rnk3,3 $F.2.txt >kk3
sort -k1,1 -rnk3,3 kk3 >kk4
perl script/topscore.pl kk4 1 $F.2.sort.txt >ll
perl script/readname.pl $F.2.sort.txt $FILE.head.txt $F.2.sort.name.txt
sed -i 's/'.param.txt'//g' $F.2.sort.name.txt
awk  '{gsub(".fa.","\t",$0); print;}' $F.2.sort.name.txt >kkk
awk  '{gsub(".fna_","\t",$0); print;}' kkk >$F.3.sort.name.txt
perl script/readfinal.pl $F.3.sort.name.txt taxonomy.txt $F.final.classify.txt

echo "Done labeling read sets and best match using SGM["$(date)"]"


rm *.2.sort.txt *.2.sort.name.txt *.3.sort.name.txt kk3 kk4 kkk ll seq.txt seq1.txt $FILE.head2.txt $FILE.loc.txt $FILE.pos.txt 



finish_time=$(date +%s)

echo "Total Time duration of classification using SGM: $((finish_time - start_time)) secs."
