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

			[Ww]* ) if [ -d "wgm_models" ]; then
					sh WGM.sh $FILE $OUTPUT.wgm
					
					echo ">>Generating Blast output for $FILE"
					blastn -query $FILE -db mydb.fa -outfmt 7 -out $FILE.BLAST_output.txt
					
					sh SGMBL.sh $FILE $OUTPUT.wgm.2.txt $FILE.BLAST_output.txt
				else
					echo "directory with wgm_models does not exist"
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
			
					sh SGMBL.sh $FILE $OUTPUT.sgm.2.txt $FILE.BLAST_output.txt
        			else
					echo "directory with sgm_models does not exist"
				fi;break;; 

			[Ww]* ) if [ -d "wgm_models" ]; then
					sh WGM.sh $FILE $OUTPUT.wgm
					
					sh SGMBL.sh $FILE $OUTPUT.wgm.2.txt $FILE.BLAST_output.txt
				else
					echo "directory with wgm_models does not exist"
				fi;break;;  

			* ) echo "Please answer sgm or wgm.";;
    			esac
		   done;break;;
        * ) echo "Please answer Yes or No.";;
    esac
done
finish_time=$(date +%s)

echo "Total Time duration of classification: $((finish_time - start_time)) secs."







