
ulimit -s unlimited

markov_order=$1 
confid1=$2
confid2=$3
confid3=$4

start_time=$(date +%s)
read -p "Do you want to create blast database(yes) or not(no)?" YN
########segmentaing and clustering of genomes ####################################################
segstart_time=$(date +%s)

echo ">>starting segmentation."
cc script/seg_clus1.c -lm
for file in genomes/*.fna; do ./a.out "$file" $markov_order $confid1 $confid2 $confid3;done >out
rm genomes/*.seg.txt

segfinish_time=$(date +%s)
echo "segmentation completed."
echo "Time duration for segmentation: $((segfinish_time - segstart_time)) secs."
######## creating the sequence file for each cluster.############

for file in genomes/*.fna; do 
perl script/seqextract.pl $file $file.clustering2.txt tempfile.txt >LL

perl script/spilit.pl tempfile.txt;done

###############Making the probabilistic models####################################################
modstart_time=$(date +%s)
if [ -d "sgm_models" ]; then
	echo "directory sgm_models exist"
else
	echo "created directory sgm_models"
	mkdir sgm_models
fi

echo ">>starting making model."
cc script/model.c -lm
for file in genomes/*.fa; do ./a.out "$file";done >out
mv genomes/*.param.txt sgm_models 

modfinish_time=$(date +%s)
echo "model generated."
echo "Time duration for model generation: $((modfinish_time - modstart_time)) secs."

rm genomes/*clustering2.txt
rm genomes/*fa
rm tempfile.txt LL seq.txt seq1.txt

##################################################################################################
##########################Making Blast database ##################################################
blaststart_time=$(date +%s)
case $YN in
	[Yy]* )	echo ">> Generating blastdatabase"
		cat genomes/*fna > mydb.fa
		makeblastdb -in mydb.fa -dbtype nucl
		echo ">>blastdatabase generated"
		break;;
	[Nn]* ) exit;;
	* ) echo "Please answer yes or no.";;
    esac
blastfinish_time=$(date +%s)
echo "Time duration for making blast database: $((blastfinish_time - blaststart_time)) secs."
####################################################################################################
finish_time=$(date +%s)
echo "Time duration for model generation and blast database: $((finish_time - start_time)) secs."
####################################################################################################






