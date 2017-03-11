#sh automate.sh test.fa output
FILE=$1
OUTPUT=$2

start_time=$(date +%s)
while true; do
    read -p "Do you want to blast testfile again(yes) or use existing(no)?" YN
    case $YN in
	[Yy]* )
		while true; do
    		read -p "Do you wish to use SGM or WGM?" SW
    		case $SW in
        		[Ss]* ) if [ -d "sgm_models" ]; then
					sh SGM.sh $FILE $OUTPUT.sgm
					
					echo ">>Generating Blast output for $FILE"
					blastn -query $FILE -db mydb.fa -outfmt 7 -out $FILE.BLAST_output.txt
					sh SGMBL.sh $FILE $OUTPUT.sgm.2.txt $FILE.BLAST_output.txt
        				
					else
						echo "directory with sgm_models does not exist"
				fi;break;;

			[Ww]* ) echo "Please enter the MARKOV MODEL ORDER:"
				read ORDER
					
				if [ -d "wgm_models" ]; then
					sh WGM.sh $FILE $OUTPUT.wgm.$ORDER $ORDER
					
					echo ">>Generating Blast output for $FILE"
					blastn -query $FILE -db mydb.fa -outfmt 7 -out $FILE.BLAST_output.txt
					
					if [ -f "$OUTPUT.wgm.$ORDER.2.txt" ]; then
						sh SGMBL.sh $FILE $OUTPUT.wgm.$ORDER.2.txt $FILE.BLAST_output.txt
					else
						echo "directory wgm_models with model of order $ORDER do not exist."
						exit
					fi

				else
					echo "directory wgm_models does not exist"
					exit
				fi;break;; 
         
				* ) echo "Please answer sgm or wgm.";;
    			esac
		done;break;;
   
	[Nn]* )
		while true; do
    		read -p "Do you wish to use SGM or WGM?" SW
    		case $SW in
        		[Ss]* ) if [ -d "sgm_models" ]; then
					sh SGM.sh $FILE $OUTPUT.sgm
			
					if [ -f "$FILE.BLAST_output.txt" ]; then
						sh SGMBL.sh $FILE $OUTPUT.sgm.2.txt $FILE.BLAST_output.txt

					else 
						echo "$FILE.BLAST_output.txt do not exist"
						exit
					fi

        			else
					echo "directory with sgm_models does not exist"
				fi;
				break;; 

			[Ww]* ) echo "Please enter the MARKOV MODEL ORDER:"
				read ORDER
					if [ -d "wgm_models" ]; then
					sh WGM.sh $FILE $OUTPUT.wgm.$ORDER $ORDER
					
						if [ -f "$OUTPUT.wgm.$ORDER.2.txt" ] && [ -f "$FILE.BLAST_output.txt" ]; then
							sh SGMBL.sh $FILE $OUTPUT.wgm.$ORDER.2.txt $FILE.BLAST_output.txt
	
							elif [ ! -f "$OUTPUT.wgm.$ORDER.2.txt" ]; then
								echo "directory wgm_models with model of order $ORDER do not exist"
								exit
							else 
								echo "$FILE.BLAST_output.txt do not exist"
								exit
							
							
						fi
					else
						echo "directory wgm_models does not exist"
						exit
					fi;break;; 

				* ) echo "Please answer sgm or wgm.";;
    				esac
		   	done;break;;
        * ) echo "Please answer Yes or No.";;
    esac
done
finish_time=$(date +%s)

mkdir Results
mv *final.classify.txt Results
mv *formulabased.taxonomy.txt Results
mv *singularhit.taxonomy.txt Results
mv *BLAST_output.txt.2.txt Results
rm $OUTPUT*
echo "Total Time duration of classification: $((finish_time - start_time)) secs."







