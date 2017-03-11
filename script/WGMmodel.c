/*
PROGRAM FOR Markov Model Segmentation and Clustering of order second when oligonucleotides not being integer..
at dof 31.5 with normal clustering
usage:
	cc seg_clus.c -lm
 	for file in seq/seq_*; do ./a.out "$file" 2 0.40 0.90 0.90; done >out
g1:file after segmentation, g2=file after contiguous cluster;g3=file after non-contiguous cluster;g4=FINAL CLUSTER FILE,,that further need to arrange by ascending order.
*/

#include<stdio.h>
#include<math.h>
#include<string.h>
#include <stdlib.h> 
#include <time.h>
#include <stddef.h>


//global variables

int ttl_oligos=0,ttl_oligos1=0,hash[200000000]={},length=0,length1=0,i=0;

FILE *fsegmlist, *fconsegmlist, *fclusterlist, *fnonconlist;

int main(int argc, char *argv[])
{
	FILE *input, *Seqfile,*Seqfile1,*fp;
	
	int cindex, rindex,sum=0,n2,z,zindex;
	int i0=0,i1=0,i2=0,i3=1,i4=0,j5=0,l0=0,g2,exp=0,n1=0,n,j1=0,j2=0,j3=0,j4=0,j0=0,k0,count,size=0,temp,k;
	char ch, line[50000],opfile[50],opfile1[20],opfile2[20]; 
       
        double entr =0.0,sumoligos=0,entr_sum=0.0;
	double *freqOligos=malloc(262144*sizeof(double)),*prob_oligos=malloc(262144*sizeof(double));	
	char *genome=malloc(200000000*sizeof(char)),opf_extn[]=".param.txt",opf_extn1[]="_1.fa",opf_extn2[]=".clustering2.txt";

	input = fopen(argv[1],"r");
 	rindex=atoi(argv[2]);
/*###################################################Step I:FILE INPUT#################################################################################################*/

     Seqfile = fopen("seq.txt","w");  

     i0=0;i1=0;
	while(1)
	{
		ch=fgetc(input);
		if(ch==EOF)
			break;
		
		if ((i1==2)&&(ch!='\n'))
		{
			if ((ch=='A')||(ch=='T')||(ch=='C')||(ch=='G')||(ch=='a')||(ch=='t')||(ch=='c')||(ch=='g')||(ch=='N')||(ch=='n'))
			{
				fputc(ch,Seqfile);
				genome[i0]=ch;
				i0++;
			}
			//else {printf ("Error ! : The input sequence contains character(s) other than A/T/C/G/N. Aborting sequence upload.\n");exit(1);}
		}		
		else if (i1==1)
		{
			if (ch=='\n') {i1=2;putchar(ch);}
			else {putchar (ch);}
		}
		else if (i1==0)
		{
			if(ch=='>'){i1=1;}
			else {i1=2;printf ("The input sequence file has no FASTA header line.\n");}
		}
	}
	

	printf ("#Sequence successfully uploaded.\n#Length of input sequence:%d\n",i0);
	fclose(input);
	fclose(Seqfile);

/*#################################################################Transitional Probablities##########################################################################################*/
 	
    

	sprintf(opfile,"%s%s.%d%s",argv[1],opf_extn1,rindex,opf_extn);
	fsegmlist=fopen(opfile,"w");

	printf ("#Order of Markov Model: %d\t",rindex);
	rindex++;
	ttl_oligos=pow(4.0,rindex);
	printf ("#Count of all possible oligomers: %d\n",ttl_oligos);
	
	
/*Converts the input (symbolic) sequence into numeric array where consecutive overlapping oligomers are mapped uniquely into a number ranging between 0 to 4^(m+1). */
 
	Seqfile = fopen("seq.txt","r");

	for (i0=0;1;i0++)
	//for (i0=69860;i0<104860;i0++)
	{
		i2=0;
		for (i1=0; i1 < rindex; i1++)				//Inner loop which encodes an oligomer into a number
		{
			ch = fgetc(Seqfile); //putchar(ch);
			if((ch=='N')||(ch=='n')){i2=1;}
			else if (i2!=1)
			{
				if ((ch == 'A') || (ch == 'a')) cindex = 0;
				if ((ch == 'T') || (ch == 't')) cindex = 1;
				if ((ch == 'G') || (ch == 'g')) cindex = 2;
				if ((ch == 'C') || (ch == 'c')) cindex = 3;
				
				exp=rindex-i1-1;
				sum = sum + cindex * pow(4.0,exp);
			}
		}
		if (i2==1){sum=ttl_oligos;}
		
		hash[i0] = sum; sum = 0; length++;
		if (fgetc(Seqfile)== EOF) break;
		fseek (Seqfile, -rindex, SEEK_CUR);
	}
	printf("#Length of input sequence:%d\n",length);	  
	fclose (Seqfile);

	//strcpy (opfile,argv[1]);
  	//strcat(argv[1],opf_extn);
	



	for(k=0;k<262144;k++)//sets the values of arrays "freqOligos" & "freqOligos1" null.
	{
		freqOligos[k]=0;
	}	
		
	n=0;

	for(k=0;k<length;k++)
	{
		temp=hash[k];
		if (temp!=ttl_oligos) {freqOligos[temp]++;n++;}
	}

	
	//for(k=0;k<64;k++)//sets the values of arrays "freqOligos" & "freqOligos1" null.
	//{
		
		for (j1=0;j1<ttl_oligos;j1++)
		{
			sumoligos=(freqOligos[j3]+freqOligos[j3+1]+freqOligos[j3+2]+freqOligos[j3+3]+4);
 			//if (sumoligos>0)
 			//{
				prob_oligos[j1]=(double)(freqOligos[j1]+1)/sumoligos;
				
				prob_oligos[j1]=log (prob_oligos[j1]);
				
				fprintf (fsegmlist,"%d\t%lf\n",j1,prob_oligos[j1]);
			
 			//}
			
			
			if (((j1+1)%4 ==0)&&(j1!=0))
			{j3+=4;/*printf ("\n");*/}
	       }
/*#################################################################Initial Probablities##########################################################################################*/
if(rindex>1){
		zindex=rindex-1;
       	 	ttl_oligos1=pow(4.0,zindex);
		printf ("#Count of all possible oligomers: %d\n",ttl_oligos1);
	
	
/*Converts the input (symbolic) sequence into numeric array where consecutive overlapping oligomers are mapped uniquely into a number ranging between 0 to 4^(m+1). */
 
		Seqfile = fopen("seq.txt","r");

		for (i0=0;1;i0++)
	//for (i0=69860;i0<104860;i0++)
		{
			i2=0;
			for (i1=0; i1 < zindex; i1++)				//Inner loop which encodes an oligomer into a number
			{
				ch = fgetc(Seqfile); //putchar(ch);
				if((ch=='N')||(ch=='n')){i2=1;}
				else if (i2!=1)
				{
					if ((ch == 'A') || (ch == 'a')) cindex = 0;
					if ((ch == 'T') || (ch == 't')) cindex = 1;
					if ((ch == 'G') || (ch == 'g')) cindex = 2;
					if ((ch == 'C') || (ch == 'c')) cindex = 3;
				
					exp=zindex-i1-1;
					sum = sum + cindex * pow(4.0,exp);
				}
			}
			if (i2==1){sum=ttl_oligos1;}
		
			hash[i0] = sum; sum = 0; length1++;
			if (fgetc(Seqfile)== EOF) break;
			fseek (Seqfile, -zindex, SEEK_CUR);
		}

		printf("#Length of input sequence:%d\n",length1);	  
		fclose (Seqfile);



		for(k=0;k<262144;k++)//sets the values of arrays "freqOligos" & "freqOligos1" null.
		{
			freqOligos[k]=0;
		}	
		
		n=0;

		for(k=0;k<length1;k++)
		{
			temp=hash[k];
			if (temp!=ttl_oligos) {freqOligos[temp]++;n++;}
		}

		printf("#Length of input sequence:%d\n",n);
		for (j1=0;j1<ttl_oligos1;j1++)
		{
			//sumoligos=(freqOligos[j2]+freqOligos[j2+1]+freqOligos[j2+2]+freqOligos[j2+3]);
 			//if (n>0)
 			//{
				prob_oligos[j1]=(double)(freqOligos[j1]+1)/(n+1);
				
				prob_oligos[j1]=log (prob_oligos[j1]);
			
				j4=j1+ttl_oligos;
				fprintf (fsegmlist,"%d\t%lf\n",j4,prob_oligos[j1]);
			
 			//}
						
			
	       }

	}		

	free(genome);
	return 1;	
}

/*###########################################################################################################################################################*/

/*double entropy (double *freq_oligos,double seqLen)
{
	int j1,j2=0;
	double ent=0.0,*prob_oligos=malloc(256*sizeof(double)),sumoligos=0;	
		
	for (j1=0;j1<ttl_oligos;j1++)
	{
		sumoligos=(freq_oligos[j2]+freq_oligos[j2+1]+freq_oligos[j2+2]+freq_oligos[j2+3]);
 		if (sumoligos>0)
 		{
			prob_oligos[j1]=(double)freq_oligos[j1]/sumoligos;
			fprintf (fsegmlist,"%lf\n",prob_oligos[j1]);
			//if (prob_oligos[j1]>0.0)
				//ent+=freq_oligos[j1]*(log (prob_oligos[j1]));
 		}
		
		if (((j1+1)%4 ==0)&&(j1!=0))
			{j2+=4;/*printf ("\n");*/
	/*}
	free(prob_oligos);
	return ent;
}

/*********************************************************************/

