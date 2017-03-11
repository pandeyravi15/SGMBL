#!/bin/sh
ulimit -s unlimited

markov_order=$1
start_time=$(date +%s) 
###############Making Models and Blast database for New genomes #################################################
# model_gen function declaration.
model_gen()
{
	
	echo ">>starting making WGM for order $markov_order."

		if [ -d "wgm_models" ]; then
			echo "directory wgm_models exist"
		else
			echo "created directory wgm_models"
			mkdir wgm_models
		fi

		cc script/WGMmodel.c -lm
		for file in addgenomes/*.fna; do ./a.out "$file" $markov_order;done >out
		mv addgenomes/*.param.txt wgm_models 
		echo "WGM for new genomes have been generated and added to wgm_models."
}
####################################################################################################
while true; do
    read -p "Do you also wish to update Blast database?" YN
    case $YN in
        [Yy]* ) model_gen 
		
		echo ">> Regenerating blastdatabase"
		cat genomes/*fna > mydb.fa
		makeblastdb -in mydb.fa -dbtype nucl
		echo ">>blastdatabase regenerated"
		break;;  

	[Nn]* ) model_gen 
		break;; 
		 
         
	* ) echo "Please answer YES or NO.";;
    esac
done
rm seq.txt
####################################################################################################
finish_time=$(date +%s)
echo "Time duration for model generation and blast database: $((finish_time - start_time)) secs."
####################################################################################################



