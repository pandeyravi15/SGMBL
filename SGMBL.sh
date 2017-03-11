start_time=$(date +%s)
F1=$1
F=$2

#H=$2
B=$3
echo ">>combining Blast output with $F:"


awk '{print $2 "\t" $1}' $F1.head.txt >key.txt
awk '{print $2 "\t" $2}' $F1.head.txt >name.txt
sed '/^#/ d' $B >$B.2.txt
python3 script/fragflip.py $F key.txt
python3 script/multimethod.py name.txt $B.2.txt $F.flip

#Taking only the tophit BLAST, and combining the scores using the formula:
#echo "Formula-based:"

perl script/label.pl $F.flip.results taxonomy.txt $F.formulabased.taxonomy.txt


finish_time=$(date +%s)
echo "Reads classified using formula."
echo "Time duration: $((finish_time - start_time)) secs."

#In cases where BLAST only identified a singular best match, this is automatically chosen and BLAST is given credit in the output file
start_time=$(date +%s)
#echo "BLAST singularhit:"
python3 script/singularblast_f.py $B.2.txt $F key.txt
perl script/label.pl $F.results taxonomy.txt $F.singularhit.taxonomy.txt
echo "Reads classified using singular hit method."
finish_time=$(date +%s)
echo "Time duration: $((finish_time - start_time)) secs."
finish_time=$(date +%s)
##remove extra files
rm key.txt name.txt $F.flip $F.flip.results

