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

int ttl_oligos=0,ttl_oligos1=0,hash[9000000]={},length=0,length1=0,x=0,z=0,i=0,j=0,g1=0;


FILE *fsegmlist, *fconsegmlist, *fclusterlist, *fnonconlist;

int main(int argc, char *argv[])
{
	FILE *input, *Seqfile,*Seqfile1,*fp,*fp1;
	
	int cindex, rindex,sum=0,n2,zindex,ravi,ravi2;
	int x0=0,i0=0,i1=0,i2=0,i3=1,u=0,v=0,j5=0,l0=0,g2,exp=0,n1=0,n,j1=0,j2=0,j3=0,j4=0,j0=0,k0,count,size=0,temp,k,pp1;
	int *freqOligos=malloc(262144*sizeof(double)),arr1[330000],arr3[10000],arr4[10000];
	char ch,opfile[20],opfile1[20],opfile2[20],line[500000]; 
        
        double entr =0.0,arr2[330000],score=0.0;
	double *prob_oligos=malloc(262144*sizeof(double)),sumoligos=0;	
	char *genome=malloc(9000000*sizeof(char)),opf_extn[]=".param.txt";


 	rindex=atoi(argv[3]);
        fsegmlist=fopen(argv[4],"a");
	 fp = fopen(argv[1], "r");
	fp1 = fopen(argv[5], "r");
    while (!feof(fp))  
    {
        fscanf(fp, "%d\t%lf\n", &arr1[u], &arr2[u]);
        //printf("%d\t%d\n",arr1[u],arr2[u]);
        u++;
    }        

	 while (!feof(fp1))  
    {
        fscanf(fp1, "%d\t%d\n", &arr3[v], &arr4[v]);
        //printf("%d\t%d\n",arr1[u],arr2[u]);
        v++;
    }        
/*###################################################Step I:FILE INPUT#################################################################################################*/

         input = fopen(argv[2],"r");
	
      Seqfile1 = fopen("seq1.txt","w");
       //Seqfile = fopen("seq.txt","w");

	j0=0;
	while(fgets (line, sizeof line, input ) != NULL)
	{
	
		if(line[0]!='>')
		{ 
			fputs ( line, Seqfile1 );j0++;
		
		}
	}	
	printf ("#Sequence successfully uploaded.\n#No of lines in input sequence:%d\n",j0);
	fclose(input);
	fclose(Seqfile1);
	
 	Seqfile1 = fopen("seq1.txt","r");
       	 Seqfile = fopen("seq.txt","w");

	i0=0;i1=0;
	while(1)
	{
		ch=fgetc(Seqfile1);
		if(ch==EOF)
		break;
		
	
		if ((ch=='A')||(ch=='T')||(ch=='C')||(ch=='G')||(ch=='a')||(ch=='t')||(ch=='c')||(ch=='g')||(ch=='N')||(ch=='n'))
		{
			fputc(ch,Seqfile);
			genome[i0]=ch;
			i0++;
		}
		//else {printf ("Error ! : The input sequence contains character(s) other than A/T/C/G/N. Aborting sequence upload.\n");exit(1);}
	}

	

	printf ("#Sequence successfully uploaded.\n#Length of input sequence:%d\n",i0);
	fclose(Seqfile1);
	fclose(Seqfile);
	
/*###########################################################Step:2#########################################################################################*/	
        
	printf ("#Order of Markov Model: %d\t",rindex);
	rindex++;
	ttl_oligos=pow(4.0,rindex);
	printf ("#Count of all possible oligomers: %d\n",ttl_oligos);
	
/*Converts the input (symbolic) sequence into numeric array where consecutive overlapping oligomers are mapped uniquely into a number ranging between 0 to 4^(m+1). */
        Seqfile = fopen("seq.txt","r");
	for (i0=0;1;i0++)
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

	
for(pp1=0;pp1<v;pp1++)
{
 ravi2++;score=0.0;

	for(k=0;k<262144;k++)//sets the values of arrays "freqOligos" & "freqOligos1" null.
	{
		freqOligos[k]=0;
	}	
		
	n=0;
	
	for(k=arr3[pp1];k<=arr4[pp1]-rindex;k++)
	{
		temp=hash[k];j5=hash[arr3[pp1]];
		if (temp!=ttl_oligos) {freqOligos[temp]++;n++;}
	}
	ravi=(j5/4)+ttl_oligos;
	score+=arr2[ravi];
	printf("#Length of input sequence:%d\t%d\t%lf\n",ravi,j5,arr2[ravi]);
	
	
		
		for (j1=0;j1<ttl_oligos;j1++)
		{
			
				score+=freqOligos[j1]*(arr2[j1]);
	
	       }
		
/*###########################################################################################################################################################*/

 fprintf (fsegmlist,"%d\t%s\t%lf\n",ravi2,argv[1],score);

}
/*********************************************************************/
free(genome);
	return 1;	
}
