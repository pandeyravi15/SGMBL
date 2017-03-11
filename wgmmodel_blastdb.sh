
ulimit -s unlimited

markov_order=$1
start_time=$(date +%s) 
read -p "Do you want to create blast database(yes) or not(no)?" YN
###############Making Models for whole genome #################################################
modstart_time=$(date +%s)
echo ">>starting making model for order $markov_order."

if [ -d "wgm_models" ]; then
	echo "directory wgm_models exist"
else
	echo "created directory wgm_models"
	mkdir wgm_models
fi

cc script/WGMmodel.c -lm
for file in genomes/*.fna; do ./a.out "$file" $markov_order;done >out
mv genomes/*.param.txt wgm_models 

rm seq.txt
modfinish_time=$(date +%s)
echo "model generated."
echo "Time duration for model generation: $((modfinish_time - modstart_time)) secs."


##################################################################################################
###############Making Blast database #############################################################
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
