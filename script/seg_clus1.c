/*
PROGRAM FOR Markov Model Segmentation and Clustering of order second when oligonucleotides not being integer..
at dof 31.5 with normal clustering
usage:
	cc seg_clus.c -lm
 	for file in seq/temp*; do ./a.out "$file" g1 g2 g3 g4; done >out
g1:file after segmentation, g2=file after contiguous cluster;g3=file after non-contiguous cluster;g4=FINAL CLUSTER FILE,,that further need to arrange by ascending order.
*/

#include<stdio.h>
#include<math.h>
#include<string.h>
#include <stdlib.h> 
#include <time.h>
#include <stddef.h>


#define ITMAX 9000
#define EPS 3.0e-7
#define FPMIN 1.0e-30 

//function declaration
double entropy (double *,double);
float segment (int *,int,int);
float gammp(float , float );
float gammln(float xx);
float cluster (int *,int,int *);
float Noncluster (int *,int *, int *, double *,int,int,int,int,int);
float Noncontiguous(int *,int *,double(*)[10000],int,int *,int *,int *,int);
//global variables
int ttl_oligos=0,hash[50000000]={},length=0,x=0,z=0,v=0,h=0,h3=0,i=0,j=0,g1=0;
int u=0,right_end=0,z0=0,seg_points[100000000]={},fragments[30000000]={},group1[10000000]={},group2[10000000]={},r[30000000],a[300000],gp[10000000]={},gp1[10000000]={},R[3000000]={},R1[3000000]={};
double sargmax=0.0,sarg_max=0.0;
double confid,confid1,confid2,confid3;			//set the CONFIDENCE LEVEL

FILE *fsegmlist, *fconsegmlist, *fclusterlist, *fnonconlist;

int main(int argc, char *argv[])
{
	FILE *input, *Seqfile,*Seqfile1,*fp;
	
	int cindex, rindex,sum=0,sizemin=0,sizemax=0,q=1,n2,y=1,w=1,z;
	double mat[10000][10000],dstrib_oligos[256]={},dstrib_oligos1[256]={},dstrib_oligos2[256]={};
	int x0=0,i0=0,i1=0,i2=0,i3=1,i8=0,i9=0,i4=0,i5=1,i7=1,j4=1,j5=0,l0=0,g2,exp=0,n1=0,n,j1=0,j6=0,j7=0,j0=0,k0,count,size=0;
	int group[10000000]={};
	char ch, line[5000],opfile[20],opfile1[20],opfile2[20]; 
        
        double enta=0.0,entb=0.0,entab=0.0,weight1=0.0,weight2=0.0;
	double jsdiv=0.0,chi_stat=0.0,dof=0.0,percent=0.0;
	double i6=0.0,entr_sum=0.0;	
	char *genome=malloc(50000000*sizeof(char)),opf_extn[]=".seg.txt",opf_extn1[]=".clustering1.txt",opf_extn2[]=".clustering2.txt";

	//rindex=2;

/* OUTPUT FILES*/
 	rindex=atoi(argv[2]);
       //fsegmlist=fopen(argv[3],"w");
	//fclusterlist=fopen(argv[4],"w");
	//fnonconlist = fopen(argv[4],"a");

/* OConfidence threshold*/	
	sscanf(argv[3],"%lf\n",&confid) ;				// for segmentation
	sscanf(argv[4],"%lf\n",&confid1);				// for contiguous clustering
	sscanf(argv[5],"%lf\n",&confid2);				// for non-contiguous clustering

/*###################################################Step I:FILE INPUT#################################################################################################*/

         input = fopen(argv[1],"r");
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


	strcpy (opfile,argv[1]);
	strcat (opfile,opf_extn);
	fsegmlist=fopen(opfile,"w");
	
	/*strcpy (opfile1,argv[1]);
	strcat (opfile1,opf_extn1);
	fconsegmlist=fopen(opfile1,"w");*/

	strcpy (opfile2,argv[1]);
	strcat (opfile2,opf_extn2);
	fclusterlist=fopen(opfile2,"w");

/*######################################################**********************First step:Segementation*****************************######################################################*/
					
	segment (hash,0,length-1);	

	if(h<4)
	{h=4;
		for (i0=h-4;i0>=0;i0-=2)
 		{	i4++;
			group[i0]=1;
			fprintf (fsegmlist,"%d\t%d\t%d\t%d\n",i3,fragments[i0+1],fragments[i0],fragments[i0]-fragments[i0+1]);i3++;
		
		}
	}
	else
	{
		for (i0=h-4;i0>=0;i0-=2)
 		{
			i4++;	
			group[i0]=1;
        	}

		for (i0=h-4;i0>=0;i0-=2)
 		{
			//fprintf (fsegmlist,"%d\t%d\t%d\t%d\t%d\n",i3,fragments[i0+1],fragments[i0],fragments[i0]-fragments[i0+1],group[i0]);  //file after segmentation
			fprintf (fsegmlist,"%d\t%d\t%d\t%d\n",i3,fragments[i0+1],fragments[i0],fragments[i0]-fragments[i0+1]);	
			i3++;	
		}
	}
	
	//printf ("%d\t%d\n",h3,h);
/***********************************************************Second step:Contiguous clustering****************************************************************************************/
	cluster (hash,0,group);

	for (i0=h-4;i0>=h3;i0-=2)
 	{
		i8++;		
		group1[i0]=i8; 
		//group[i0]=1;
		size =(fragments[i0]-fragments[i0+1]);
		percent = (double)size/length;
		//fprintf (fconsegmlist,"%d\t%d\t%d\t%d\t%d\t%f\t%d\n",i8,1,fragments[i0+1],fragments[i0],fragments[i0]-fragments[i0+1],percent,group[i0]); // output cluster of contiguous clustering//
        
        }

//printf ("%d\t%d\n",h3,h);
	
/*************************************************************Third step:contiguous+NonContiguous clustering**************************************************************************/	
	for (i0=h-4;i0>=h3;i0-=2)
	{
		j1++;
		for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos1[i1]=0;
		}

		for (i1=fragments[i0+1];i1<=fragments[i0];i1++)		
 		{
			i2=hash[i1];
			dstrib_oligos1[i2]++;
			
		}
		
                for (i1=0;i1<ttl_oligos;i1++)				
 		{
			
			dstrib_oligos1[i1]=dstrib_oligos1[i1]/(group[i0]);
		}
		n1=0;
		for (i1=0;i1<ttl_oligos;i1++)				
 		{
			n1+=dstrib_oligos1[i1];
		}

			size = fragments[i0]-fragments[i0+1];
			percent = (double) size/length;
		//fprintf (fclusterlist,"\n%d\t%d\t%d\t%d\t%d\t%d\t%f\n",j1,group1[i0],fragments[i0+1],fragments[i0],fragments[i0]-fragments[i0+1],size,percent);
		
		gp[z0]=group1[i0];gp1[z0]=group[i0];z0=z0+2;	
		r[i]=fragments[i0+1];i++;
		r[i]=fragments[i0];i++;

 		Noncluster(hash,group,group1,dstrib_oligos1,j1,i0,h3,size,length);

		for (x0=j;x0<i;x0+=2)
		{
			a[x0]=j1;
		}
		j=i;

		n1=0;
		for (i1=0;i1<ttl_oligos;i1++)				
 		{
			mat[i1][i0] = dstrib_oligos1[i1];
			n1+=mat[i1][i0];
		}	
	}

	size=0;
	for (i0=0;i0<i;i0+=2)
	{	
		if(a[i0] != a[i0-2]){size=0;}
		size =size + (r[i0+1]-r[i0]);
		percent = (double)size/length;
		fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%d\t%d\t%lf\n",a[i0],r[i0],r[i0+1],gp1[i0],r[i0+1]-r[i0],size,percent);  // output cluster of non-contiguous clustering//
		//fprintf (fclusterlist,"%s\t%d\t%d\t%d\t%d\n",argv[1],a[i0],r[i0],r[i0+1],gp1[i0]);
		R[i0]=a[i0];R1[i0]=gp[i0];
	}
	
	for (i0=h-4;i0>=h3;i0-=2)
	{
		j7++;
		gp[i0]=j7;
		
		//printf ("%s\t%d\t%d\t%d\t%d\t%d\n",argv[1],gp[i0],group1[i0],fragments[i0+1],fragments[i0],group[i0]);
	}

/*#####################################**********************Final step:ReNonContiguous clustering**************************#######################################################*/
	//Noncontiguous (group,group1,mat,h3,gp,gp1,R,i);


	printf ("%d\t%d\n",h3,h);
	
	for (i0=h-4;i0>=h3;i0-=2)
	{
		//j7++;
		//gp[i0]=j7;
		
		//printf ("%s\t%d\t%d\t%d\t%d\t%d\n",argv[1],gp[i0],group1[i0],fragments[i0+1],fragments[i0],group[i0]);
	}

	for (i0=0;i0<i;i0+=2)
	{	
		if(R[i0] != R[i0-2]){size=0;}
		size =size + (r[i0+1]-r[i0]);
		percent = (double)size/length;
		//fprintf (fclusterlist,"%s\t%d\t%d\t%d\t%d\t%lf\n",argv[1],R[i0],r[i0],r[i0+1],gp1[i0],percent);  // FINAL output cluster of clustering//
	}

	free(genome);
	return 1;	

}	
//END OF main()
	
/*################################################STEP:3:: function for iterative segmentation######################################################*/

float segment (int *hash,int k1,int k2)
{
	
	int n=0,n1=0,n2=0,n3=0,pos=0,k=0,k4=0,temp,block_index,temp_var1=0,interval=0;	
	double max=0.0,entr=0.0,entr1=0.0,entr2=0.0,entr3=0.0,entr4=0.0,pr_entr1=0.0,pr_entr2=0.0,div=0.0,w1=0.0,w2=0.0,a,b,c,d,inter;
	double neff=0.0,sarg=0.0,smax=0.0,sx=0.0,beta=0.0,dof=0.0;	
	double *freqOligos=malloc(64*sizeof(double)),*freqOligos1=malloc(64*sizeof(double));
	
	
	//printf ("Input string: start:%d end:%d size=%d\n",k1,k2,k2-k1);
	
	for(k=0;k<64;k++)//sets the values of arrays "freqOligos" & "freqOligos1" null.
	{
		freqOligos[k]=0;freqOligos1[k]=0;
	}	
		
	n=0;
	for(k=k1;k<=k2;k++)
	{
		temp=hash[k];
		if (temp!=ttl_oligos) {freqOligos[temp]++;n++;}
	}

	inter =(n/10000);
	if (inter > 1){interval= inter;}
	else{interval = 1 ;}
	
	if (n>1024)
	{
		
		entr=entropy(freqOligos,n);		
		entr=  (1.0/n) * entr ;
		entr= (-entr)/log (2.0);
			
			  
		for(k=k1;k<=k2-15;k++)
		{	
			temp=hash[k];			
		    if ((temp!=ttl_oligos) && (n1<n-1) )
		    { 
			
			freqOligos1[temp]++; 
			freqOligos[temp]--;
												
			n1++;n2=n-n1;	
	
			//if((n1>=256)&&(n1%(500)==0))
			if((n1>=15)&&(n1%(interval)==0))
			{			
				entr1=entropy(freqOligos1,n1);			
				entr1=  (1.0/n1) * entr1 ;
				entr1= (-entr1)/log (2.0);	
				
						
				entr2=entropy(freqOligos,n2);			
				entr2=  (1.0/n2) * entr2 ;
				entr2= (-entr2)/log (2.0);	


				w1=(double)n1/n; w2=(double)n2/n;
				div=entr-(w1*entr1)-(w2*entr2);
			
				if (div>max){max=div;pos=k;entr3=entr1;entr4=entr2;}
			
			}
 			
		    }
		}
	
	
if (ttl_oligos==4)	
		{  a=2.7784; b=-7.97084; c=0.0; d=0.80;dof=1.5;	}		
                        						
else if (ttl_oligos==16)
		{  a=1.55681; b=-2.19548; c=0.0; d=0.946; dof=7.5;}
							
else if (ttl_oligos==64)
		{  a=1.13049; b=-2.44732;c=0.002304; d=1.025173;dof=31.5;}
		// {  a=2.39; b=-7.66;c=0.0029; d=0.841;}
				
                           beta=(c*log(n)) + d;
                           neff=(a*log(n)) + b;
		           sarg=(log(2.0)*beta*n*max);
			   sx=gammp(dof,sarg);
			   //sx=gammp(24.0,sarg);
			   smax=pow(sx,neff);
			
//condition to check the significance of D_max
		if (smax>confid)
		{	
			//printf ("%lf\n",sarg);
			temp_var1=pos;
			//condition where sequence is segmented for first time
			if(u==0)
			{	
				u++;
				seg_points[v]=k1;
				v++;		
				seg_points[v]=temp_var1;
				v++;	
				seg_points[v]=temp_var1;
				v++;	
				seg_points[v]=k2;
				
				//printf("**Segment [%d...%d] partitioned at %d with  signif: %lf\n",k1,k2,temp_var1,smax);
				
				segment(hash,seg_points[v-1],seg_points[v]);
				
			}
//condition where segments are segmented subsequently
			else
			{	
	
					seg_points[v-1]=k1;
					seg_points[v]=temp_var1;
					v++;
					seg_points[v]=temp_var1;
					v++;
					seg_points[v]=k2;
					
					//printf("**Segment [%d...%d] partitioned at %d with  signif: %lf\n",k1,k2,temp_var1,smax);
				
                                       segment(hash,seg_points[v-1],seg_points[v]);

			}
		}
		//condition when the D_max is not significant 
		else
		{
			
			//printf("<< Not significant>> Loc: %d...%d  at <%d> Signf: %lf size=%d JSdiv_max=%lf\n",seg_points[v-1],seg_points[v],pos,smax,seg_points[v]-seg_points[v-1],max);

			if(v>0)
			{
 				fragments[h]=seg_points[v];
				h++;
				fragments[h]=seg_points[v-1];
				h++;
				right_end=seg_points[v];	
			
				out:
				v=v-2;
				//printf ("\nfun called for %d--%d #fragments=%d\n",seg_points[v-1],seg_points[v],h);
                        	segment(hash,seg_points[v-1],seg_points[v]);
			}

			else
			{

				fragments[h]=k2;
				h++;
				fragments[h]=k1;
				h++;
				right_end=seg_points[v];	
				//printf("<< Not significant>> Loc: %d...%d  at <%d> Signf: %lf size=%d JSdiv_max=%lf\n",k1,k2,pos,smax,k2-k1,max);

			}

		}
	}
	else 
	{
		free(freqOligos);free(freqOligos1);	
		
		fragments[h]=k2;
		h++;
		fragments[h]=k1;
		h++;
		right_end=seg_points[v];
		
		//printf("<< Size below threshold >> Loc: %d...%d size::%d #fragments=%d\n",k1,k2,(k2-k1),h);
		if (v>0){goto out;}
	}
	//return;
}

/**#############################**********************function for iterative contiguous clustering************##################################################################*/
float cluster(int *hash,int h2,int *group)
{
         int i1=0,i2=0,i3=1,i5=1,i7=1,i0=1,p=1,q=1,l0=0,r=group[h-4];
	double enta=0.0,entb=0.0,entab=0.0,weight1=0.0,weight2=0.0,a, b,c, d;
	double jsdiv=0.0,chi_stat=0.0,dof=0.0,signif3=0.0,neff=0.0,sx=0.0,beta=0.0;	
	double dstrib_oligos[256]={},dstrib_oligos1[256]={},dstrib_oligos2[256]={},n1=0,n2=0,n=0;

	for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos1[i1]=0;
		}
		for (i1=fragments[h-3];i1<=fragments[h-4];i1++)		
 		{
			i2=hash[i1];
			dstrib_oligos1[i2]++;

		}
                for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos1[i1]=dstrib_oligos1[i1]/group[h-4];
		}
		n1=0;
		for (i1=0;i1<ttl_oligos;i1++)				
 		{
			n1+=dstrib_oligos1[i1];
		}
	p=group[h-4];
	//fprintf (fsegmlist,"%d\n",n1);
	h2=h3;
	for (i0=h-6;i0>=h2;i0-=2)
 	{i5++;
                for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos2[i1]=0;
		}
		for (i1=fragments[i0+1];i1<=fragments[i0];i1++)		
 		{
			i2=hash[i1];
			dstrib_oligos2[i2]++;
		
                 }
                for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos2[i1]=dstrib_oligos2[i1]/group[i0];
		}

                n2=0;
                for (i1=0;i1<ttl_oligos;i1++)				
 		{
			n2+=dstrib_oligos2[i1];
		}
if(n1>0){
		enta=entropy (dstrib_oligos1,n1);
		enta=(1.0/(n1)) * enta ;
		enta=(-enta)/log(2.0);

		
		entb=entropy (dstrib_oligos2,n2);
		entb=(1.0/(n2)) * entb ;
		entb=(-entb)/log(2.0);
			
		for (i1=0;i1<ttl_oligos;i1++)
                {
                  dstrib_oligos[i1]=0;
                }

		for (i1=0;i1<ttl_oligos;i1++)
                {
		dstrib_oligos[i1]=dstrib_oligos1[i1]+dstrib_oligos2[i1];                 
		}

                n=n1+n2;
		entab=entropy (dstrib_oligos,n);
		entab=  (1.0/(n)) * entab ;
		entab= (-entab)/log (2.0);

                weight1=(double)(n1)/(n);
		weight2=(double)(n2)/(n);

		
		jsdiv = entab - weight1*enta - weight2*entb ;

if (ttl_oligos==4)	
		{  a=2.7784; b=-7.97084; c=0.0; d=0.80;dof=1.5;	}		
                        						
else if (ttl_oligos==16)
		{  a=1.55681; b=-2.19548; c=0.0; d=0.946; dof=7.5;}
							
else if (ttl_oligos==64)
		{  a=1.13049; b=-2.44732;c=0.002304; d=1.025173;dof=31.5;}

		beta=(c*log(n)) + d;
                neff=(a*log(n)) + b;		  
		chi_stat = log(2.0)*(n)*jsdiv*beta;
		sx = gammp(dof,chi_stat);
		//sx = gammp(24.0,chi_stat);
		signif3 =pow(sx,neff);
//#####################################################################################################################################################/
if(signif3 < confid1)
  {
         //fprintf (fconsegmlist,"%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%d\n",i3,i5,signif3,n1,n2,n,dstrib_oligos[1],group[i0],p);
         	i7++;
      	 	for (i1=0;i1<ttl_oligos;i1++)
     	 	{
       	  	    dstrib_oligos1[i1]=((dstrib_oligos1[i1]*p)+(dstrib_oligos2[i1]*group[i0]))/(p+group[i0]);
			//dstrib_oligos1[i1]=(dstrib_oligos1[i1]+dstrib_oligos2[i1]);
         	}

     		p=p+group[i0];
     		q=q+1;
     		n1=0;

     		for (i1=0;i1<ttl_oligos;i1++)				
    		{
	  	  n1=n1+dstrib_oligos1[i1];
     		}

     		group[i0]=p;
  }

else{
		p=group[i0];

		//fprintf (fconsegmlist,"%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%d\n",i3,i5,signif3,n1,n2,n,dstrib_oligos[1],group[i0],p);
        	n1=0;
        	for (i1=0;i1<ttl_oligos;i1++)
        	{
          	  dstrib_oligos1[i1]=dstrib_oligos2[i1];
        	}
        	for (i1=0;i1<ttl_oligos;i1++)				
        	{
	  	   n1=n1+dstrib_oligos1[i1];
        	} 	  
    }
}

else{
        	p=group[i0];
		//fprintf (fconsegmlist,"%d\t%d\t%f\t%f\t%f\t%f\t%f\t%d\t%d\n",i3,i5,signif3,n1,n2,n,dstrib_oligos[1],group[i0],p);
       		 
        	for (i1=0;i1<ttl_oligos;i1++)
        	{
          	  dstrib_oligos1[i1]=dstrib_oligos2[i1];
        	}

		n1=0;
        	for (i1=0;i1<ttl_oligos;i1++)				
        	{
	  	  n1=n1+dstrib_oligos1[i1];
                }
     }

         	fragments[i0+2*(q-1)]=fragments[i0];
         	fragments[i0+2*(q-1)-1]=fragments[i0-1];
             	group[i0+2*(q-1)]=group[i0];

         	i3++;
}
		h3=h3+2*(i7-1); 

//fprintf (fsegmlist,"%d\t%d\n",i7,h3);

  		if(i7>1)
  		{
   		  cluster(hash,h3,group);
  
                } 
}

/*##################################################**Function for Contiguous & Non-contiguous Clustering***************#######################################*/

float Noncluster(int *hash,int *group,int *group1,double *dstrib_oligos1,int j1,int i0,int h4,int size,int length)
{
         int   i1=0,i2=0,l0=0,j0=0,j4=1,s=1,z1=0,w=1,w2=0,w1=0,z2=0;
	double enta=0.0,entb=0.0,entab=0.0,weight1=0.0,weight2=0.0,percent = 0.0,neff=0.0,sx=0.0,beta=0.0,a, b,c, d;
	double jsdiv=0.0,chi_stat=0.0,dof=0.0,signif4=0.0,signif6=0.0;	
	double dstrib_oligos[256]={},dstrib_oligos2[256]={},n1=0,n2=0,n;


		n1=0;
		for (i1=0;i1<ttl_oligos;i1++)				
 		{
			n1+=dstrib_oligos1[i1];
		}
//fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%f\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,group[i0]);
//fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%d\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,group[i0]);

h4=h3;
s=1;
z1=0;
	for (l0=i0-2;l0>=h4;l0-=2)
	{

		for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos2[i1]=0;
		}
		for (i1=fragments[l0+1];i1<=fragments[l0];i1++)		
 		{
			i2=hash[i1];
			dstrib_oligos2[i2]++;
		
                 }
                for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos2[i1]=(dstrib_oligos2[i1])/(group[l0]);
		}

                n2=0;
                for (i1=0;i1<ttl_oligos;i1++)				
 		{
			n2+=dstrib_oligos2[i1];
		}

if((n1>0)&&(n2>0))
{
		enta=entropy (dstrib_oligos1,n1);
		enta=(1.0/(n1)) * enta ;
		enta=(-enta)/log(2.0);
		
		entb=entropy (dstrib_oligos2,n2);
		entb=(1.0/(n2)) * entb ;
		entb=(-entb)/log(2.0);
			
		for (i1=0;i1<ttl_oligos;i1++)
                {
                  dstrib_oligos[i1]=0;
                }

		for (i1=0;i1<ttl_oligos;i1++)
                {
		dstrib_oligos[i1]=dstrib_oligos1[i1]+dstrib_oligos2[i1];                 
		}

		n=0;
                n=n1+n2;
		entab=entropy (dstrib_oligos,n);
		entab=  (1.0/(n)) * entab ;
		entab= (-entab)/log (2.0);

                weight1=(double)(n1)/(n);
		weight2=(double)(n2)/(n);
		
		jsdiv = entab - weight1*enta - weight2*entb ;

		if (ttl_oligos==4)	
		{  a=2.7784; b=-7.97084; c=0.0; d=0.80;dof=1.5;	}		
                        						
else if (ttl_oligos==16)
		{  a=1.55681; b=-2.19548; c=0.0; d=0.946; dof=7.5;}
							
else if (ttl_oligos==64)
		{  a=1.13049; b=-2.44732;c=0.002304; d=1.025173;dof=31.5;}
		beta=(c*log(n)) + d;
                neff=(a*log(n)) + b;		  
		chi_stat = log(2.0)*(n)*jsdiv*beta;
		sx = gammp(dof,chi_stat);
		//sx = gammp(24.0,chi_stat);
		signif4 =pow(sx,neff);
/************************************************************************************************************************************************************************************************/
if(signif4 < confid2)
{		
		s=s+1;
		//fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%lf\t%d\t%d\t%d\n",j1,group1[l0],fragments[l0+1],fragments[l0],n1,group[l0],s,z1);
//fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%lf\t%f\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,signif4,group1[l0],fragments[l0+1],fragments[l0],n2,n,group[i0],s,z1);
		j4++;

     		for (i1=0;i1<ttl_oligos;i1++)
     		{
         		dstrib_oligos1[i1]=((dstrib_oligos1[i1]*group[i0])+(dstrib_oligos2[i1]*group[l0]))/(group[i0]+group[l0]);
     		}

       		group[i0]=group[i0]+group[l0];
        	
     		n1=0;
     		for (i1=0;i1<ttl_oligos;i1++)				
     		{
	    	n1=n1+dstrib_oligos1[i1];
     		}

		size =size + (fragments[l0]-fragments[l0+1]);
		percent = (double)size/length;
		//a fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%d\t%d\t%f\n",j1,group1[l0],fragments[l0+1],fragments[l0],fragments[l0]-fragments[l0+1],size,percent);
		gp[z0]=group1[l0];gp1[z0]=group[i0];z0=z0+2;	
		r[i]=fragments[l0+1];i++;
		r[i]=fragments[l0];i++;
/*##############################################################**Start:recursive clutseing by first way***####################################################################################***/
    	w=1;
      w1=l0+2*(s-1);	
      //w2=i0+2*(j1-1);

	for(j0=i0-2;j0>=w1;j0-=2)
    	{  
	  //if(j0!=i0){	
	  	for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos2[i1]=0;
		}
		for (i1=fragments[j0+1];i1<=fragments[j0];i1++)		
 		{
			i2=hash[i1];
			dstrib_oligos2[i2]++;
		
                 }
                for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos2[i1]=(dstrib_oligos2[i1])/(group[j0]);
		}

                n2=0;
                for (i1=0;i1<ttl_oligos;i1++)				
 		{
			n2+=dstrib_oligos2[i1];
		}

	if((n1>0)&&(n2>0))
	{
		enta=entropy (dstrib_oligos1,n1);
		enta=(1.0/(n1)) * enta ;
		enta=(-enta)/log(2.0);
		
		entb=entropy (dstrib_oligos2,n2);
		entb=(1.0/(n2)) * entb ;
		entb=(-entb)/log(2.0);
			
		for (i1=0;i1<ttl_oligos;i1++)
                {
                  dstrib_oligos[i1]=0;
                }

		for (i1=0;i1<ttl_oligos;i1++)
                {
		dstrib_oligos[i1]=dstrib_oligos1[i1]+dstrib_oligos2[i1];                 
		}

		n=0;
                n=n1+n2;
		entab=entropy (dstrib_oligos,n);
		entab=  (1.0/(n)) * entab ;
		entab= (-entab)/log (2.0);

                weight1=(double)(n1)/(n);
		weight2=(double)(n2)/(n);
		
		jsdiv = entab - weight1*enta - weight2*entb ;

		if (ttl_oligos==4)	
		{  a=2.7784; b=-7.97084; c=0.0; d=0.80;dof=1.5;	}		
                        						
else if (ttl_oligos==16)
		{  a=1.55681; b=-2.19548; c=0.0; d=0.946; dof=7.5;}
							
else if (ttl_oligos==64)
		{  a=1.13049; b=-2.44732;c=0.002304; d=1.025173;dof=31.5;}
		beta=(c*log(n)) + d;
                neff=(a*log(n)) + b;		  
		chi_stat = log(2.0)*(n)*jsdiv*beta;
		sx = gammp(dof,chi_stat);
		//sx = gammp(24.0,chi_stat);
		signif6 =pow(sx,neff);
/************************************************************************************************************************************************************************************************/
	if(signif6 < confid2)
	{	s=s+1;
		w=w+1;
		//z1=z1-1;

		//fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%lf\t%d\t%d\t%d\n",j1,group1[j0],fragments[j0+1],fragments[j0],n1,group[j0],s,z1);
//fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%lf\t%f\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,signif6 ,group1[j0],fragments[j0+1],fragments[j0],n2,n,group[j0],s,z1);
		j4++;

     		for (i1=0;i1<ttl_oligos;i1++)
     		{
         		dstrib_oligos1[i1]=((dstrib_oligos1[i1]*group[i0])+(dstrib_oligos2[i1]*group[j0]))/(group[i0]+group[j0]);
     		}

       		group[i0]=group[i0]+group[j0];
        	
     		n1=0;
     		for (i1=0;i1<ttl_oligos;i1++)				
     		{
	    		n1=n1+dstrib_oligos1[i1];
     		}

		size =size + (fragments[j0]-fragments[j0+1]);
		percent = (double)size/length;
		//a fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%d\t%d\t%f\n",j1,group1[j0],fragments[j0+1],fragments[j0],fragments[j0]-fragments[j0+1],size,percent);
		gp[z0]=group1[j0];gp1[z0]=group[i0];z0=z0+2;
		r[i]=fragments[j0+1];i++;
		r[i]=fragments[j0];i++;
	}
/************************************************************************************************************************************************************************************************/
	else
           {z2++;
//fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%lf\t%f\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,signif6,group1[j0],fragments[j0+1],fragments[j0],n2,n,group[i0],s,z1);

		group1[j0+2*(w-1)]=group1[j0];
	 	fragments[j0+2*(w-1)]=fragments[j0];
       		fragments[j0+2*(w-1)+1]=fragments[j0+1];
             	group[j0+2*(w-1)]=group[j0];
	    }
        }
    
        else
           {z2++;      
//fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%lf\t%f\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,signif6,group1[j0],fragments[j0+1],fragments[j0],n2,n,group[i0],s,z1);
		group1[j0+2*(w-1)]=group1[j0];
	 	fragments[j0+2*(w-1)]=fragments[j0];
       		fragments[j0+2*(w-1)+1]=fragments[j0+1];
             	group[j0+2*(w-1)]=group[j0];
	   }
         }
    //  }

/*##############################################################**Stop:recursive clutseing by 1st way***####################################################################################***/
}
else
    {
		z1++;
//fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%lf\t%f\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,signif4,group1[l0],fragments[l0+1],fragments[l0],n2,n,group[i0],s,z1); 

	    	group1[l0+2*(s-1)]=group1[l0];
	 	fragments[l0+2*(s-1)]=fragments[l0];
       		fragments[l0+2*(s-1)+1]=fragments[l0+1];
             	group[l0+2*(s-1)]=group[l0];
        
    }

		//fprintf (fsegmlist,"%d\t%d\t%f\t%d\t%d\t%d\t%d\n",j1,j2,signif3,n1,n2,n,dstrib_oligos[1]);
}
else
   {
       		z1++;
//fprintf (fclusterlist,"%d\t%d\t%d\t%d\t%lf\t%f\t%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,signif4,group1[l0],fragments[l0+1],fragments[l0],n2,n,group[i0],s,z1); 

		group1[l0+2*(s-1)]=group1[l0];
	 	fragments[l0+2*(s-1)]=fragments[l0];
       		fragments[l0+2*(s-1)+1]=fragments[l0+1];
             	group[l0+2*(s-1)]=group[l0];
       
   }
}

	h3=h3+2*(s-1);

           		/*######**Start:recursive clutseing by 2nd way***#########################################***/
	       		/*if(s>1)
	       		{	
	        		Noncluster(hash,group,group1,dstrib_oligos1,j1,i0,h3,size,length);
	       		}

            		/********Stop:recursive clutseing by 2nd way***#########################################**/

if(percent > 0.60){printf ("%d\t%f\n",j1,percent);}
//else {printf ("No cluster contains atleast 60 percent of native genome \n");}

}  
/*#######################################**Stop:recursive clutseing ***####################################################################################***/

float Noncontiguous (int *group,int *group1,double (*mat)[10000],int h4,int *gp,int *gp1,int *R,int i)
{
         int i1=0,l0=0,n5=0,y=1,z=0,j0=0,s=1,z1=0,w=1,w1=0,z2=0,i0=0,j5=0,F=0,x=0;
	double enta=0.0,entb=0.0,entab=0.0,weight1=0.0,weight2=0.0,a,b,c,d;
	double jsdiv=0.0,chi_stat=0.0,dof=0.0,signif5=0.0,signif7=0.0,sx=0.0,beta=0.0,neff=0.0;	
	double dstrib_oligos[256]={},dstrib_oligos1[256]={},dstrib_oligos2[256]={},n1=0,n2=0,n=0;
	
for (i0=h-4;i0>=h3;i0-=2)
{
		j5++;

		for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos1[i1]=0;
		}
	
		for (i1=0;i1<ttl_oligos;i1++)				
 		{
			dstrib_oligos1[i1]=mat[i1][i0];			 			
		}

		n1=0;
		for (i1=0;i1<ttl_oligos;i1++)				
 		{		 
		n1+=dstrib_oligos1[i1]; 
		}
		F=0;
		//printf ("%d\t%d\t%d\t%d\t%d\t%lf\t%d\n",j5,gp[i0],group1[i0],fragments[i0+1],fragments[i0],n1,group[i0]);F=gp[i0];gp[i0]=j5;
		for (x=0;x<=i;x+=2)
		{ if(R[x]==F){R[x]=j5;gp1[x]=group[i0];}}

y=1;
for (l0=i0-2;l0>=h3;l0-=2)
{
		for (i1=0;i1<ttl_oligos;i1++)				
 		{
		  dstrib_oligos2[i1]=0;
		}

		n2=0;
		for (i1=0;i1<ttl_oligos;i1++)				
 		{
			
		  dstrib_oligos2[i1]=mat[i1][l0];
		  n2+=dstrib_oligos2[i1];
		}
//fprintf (fsegmlist,"%d\t%d\t%d\t%d\t%d\n",j5,group1[l0],fragments[l0+1],fragments[l0],group[l0]);
//}
if((n1>0)&&(n2>0))
{
		enta=entropy (dstrib_oligos1,n1);
		enta=(1.0/(n1)) * enta ;
		enta=(-enta)/log(2.0);

		
		entb=entropy (dstrib_oligos2,n2);
		entb=(1.0/(n2)) * entb ;
		entb=(-entb)/log(2.0);

		for (i1=0;i1<ttl_oligos;i1++)
                {
                  dstrib_oligos[i1]=0;
                }

		for (i1=0;i1<ttl_oligos;i1++)
                {
		dstrib_oligos[i1]=dstrib_oligos1[i1]+dstrib_oligos2[i1];                 
		}
		
		n=0;
                n=n1+n2;
		entab=entropy (dstrib_oligos,n);
		entab=  (1.0/(n)) * entab ;
		entab= (-entab)/log (2.0);

                weight1=(double)(n1)/(n);
		weight2=(double)(n2)/(n);

		
		jsdiv = entab - weight1*enta - weight2*entb ;
		if (ttl_oligos==4)	
		{  a=2.7784; b=-7.97084; c=0.0; d=0.80;dof=1.5;	}		
                        						
else if (ttl_oligos==16)
		{  a=1.55681; b=-2.19548; c=0.0; d=0.946; dof=7.5;}
							
else if (ttl_oligos==64)
		{  a=1.13049; b=-2.44732;c=0.002304; d=1.025173;dof=31.5;}
		beta=(c*log(n)) + d;
                neff=(a*log(n)) + b;		  
		chi_stat = log(2.0)*(n)*jsdiv*beta;
		sx = gammp(dof,chi_stat);
		signif5 =pow(sx,neff);
/************************************************************************************************************************************************************************************************/
if(signif5 < confid3)
{	
F=0;
		//fprintf (fnonconlist,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n",j5,group1[i0],fragments[i0+1],fragments[i0],fragments[l0+1],fragments[l0],n1,n2,group[i0],jsdiv);
		//printf ("%d\t%d\t%d\t%d\t%d\t%d\n",j5,gp[l0],group1[l0],fragments[l0+1],fragments[l0],group[l0]); F=gp[l0];
		z++;
		
		y=y+1;
		
		for (i1=0;i1<ttl_oligos;i1++)				
 		{
	         dstrib_oligos1[i1]=(dstrib_oligos1[i1]*group[i0] + dstrib_oligos2[i1]*group[l0])/(group[i0]+group[l0]) ;
			mat[i1][i0]=dstrib_oligos1[i1];
		}

		group[i0]=group[i0]+group[l0];

		n1=0;
		for (i1=0;i1<ttl_oligos;i1++)				
     		{
		n1=n1+dstrib_oligos1[i1];
     		}	 
			
		for (x=0;x<=i;x+=2)
		{ if(R[x]==F){R[x]=j5;gp1[x]=group[i0];}}
/*##############################################################**Start:recursive clutseing by first way***####################################################################################***/
      w=1;
      w1=l0+2*(y-1);
      for(j0=i0-2;j0>=w1;j0-=2)
	  {	
	  	for (i1=0;i1<ttl_oligos;i1++)				
 		{
		  dstrib_oligos2[i1]=0;
		}

		n2=0;
		for (i1=0;i1<ttl_oligos;i1++)				
 		{
			
		  dstrib_oligos2[i1]=mat[i1][j0];
		  n2+=dstrib_oligos2[i1];
		}

	if((n1>0)&&(n2>0))
	{
		enta=entropy (dstrib_oligos1,n1);
		enta=(1.0/(n1)) * enta ;
		enta=(-enta)/log(2.0);
		
		entb=entropy (dstrib_oligos2,n2);
		entb=(1.0/(n2)) * entb ;
		entb=(-entb)/log(2.0);
			
		for (i1=0;i1<ttl_oligos;i1++)
                {
                  dstrib_oligos[i1]=0;
                }

		for (i1=0;i1<ttl_oligos;i1++)
                {
		dstrib_oligos[i1]=dstrib_oligos1[i1]+dstrib_oligos2[i1];                 
		}

		n=0;
                n=n1+n2;
		entab=entropy (dstrib_oligos,n);
		entab=  (1.0/(n)) * entab ;
		entab= (-entab)/log (2.0);

                weight1=(double)(n1)/(n);
		weight2=(double)(n2)/(n);
		
		jsdiv = entab - weight1*enta - weight2*entb ;
		if (ttl_oligos==4)	
		{  a=2.7784; b=-7.97084; c=0.0; d=0.80;dof=1.5;	}		
                        						
else if (ttl_oligos==16)
		{  a=1.55681; b=-2.19548; c=0.0; d=0.946; dof=7.5;}
							
else if (ttl_oligos==64)
		{  a=1.13049; b=-2.44732;c=0.002304; d=1.025173;dof=31.5;}
		
		beta=(c*log(n)) + d;
                neff=(a*log(n)) + b;		  
		chi_stat = log(2.0)*(n)*jsdiv*beta;
		sx = gammp(dof,chi_stat);
		signif7 =pow(sx,neff);
/************************************************************************************************************************************************************************************************/
	if(signif7 < confid3)
	{	y=y+1;
		w=w+1;
		//z1=z1-1;
		F=0;
		//printf ("%d\t%d\t%d\t%d\t%d\t%d\n",j5,gp[j0],group1[j0],fragments[j0+1],fragments[j0],group[j0]); F=gp[j0];
//fprintf (fnonconlist,"%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,signif6 ,group1[j0],fragments[j0+1],fragments[j0],n2,n,group[j0],s,z1);
		
		z++;
     		for (i1=0;i1<ttl_oligos;i1++)				
 		{
	         dstrib_oligos1[i1]=(dstrib_oligos1[i1]*group[i0] + dstrib_oligos2[i1]*group[j0])/(group[i0]+group[j0]) ;
			mat[i1][i0]=dstrib_oligos1[i1];
		}

		group[i0]=group[i0]+group[j0];

		n1=0;
		for (i1=0;i1<ttl_oligos;i1++)				
     		{
		n1=n1+dstrib_oligos1[i1];
     		}	 
		for (x=0;x<=i;x+=2)
		{ if(R[x]==F){R[x]=j5;gp1[x]=group[i0];}}
	}
/************************************************************************************************************************************************************************************************/
	else
           {z2++;
//fprintf (fnonconlist,"%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,signif6,group1[j0],fragments[j0+1],fragments[j0],n2,n,group[i0],s,z1);

		group1[j0+2*(w-1)]=group1[j0];gp[j0+2*(w-1)]=gp[j0];
  		for (i1=0;i1<ttl_oligos;i1++)
  		{mat[i1][j0+2*(w-1)]=mat[i1][j0];
		}

		fragments[j0+2*(w-1)]=fragments[j0];
		fragments[j0+2*(w-1)+1]=fragments[j0+1];
		group[j0+2*(w-1)]=group[j0];
	    }
        }
    
        else
           {z2++;       /*to count unclustered segments*/
//fprintf (fnonconlist,"%d\t%d\t%d\t%d\t%d\t%f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",j1,group1[i0],fragments[i0+1],fragments[i0],n1,signif6,group1[j0],fragments[j0+1],fragments[j0],n2,n,group[i0],s,z1);
		group1[j0+2*(w-1)]=group1[j0];gp[j0+2*(w-1)]=gp[j0];
  		for (i1=0;i1<ttl_oligos;i1++)
  		{mat[i1][j0+2*(w-1)]=mat[i1][j0];
		}

		fragments[j0+2*(w-1)]=fragments[j0];
		fragments[j0+2*(w-1)+1]=fragments[j0+1];
		group[j0+2*(w-1)]=group[j0];
	   }
         }
      }

/*##############################################################**Stop:recursive clutseing by 1st way***####################################################################################***/

else
   {z1++;
		//fprintf (fnonconlist,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n",j5,group1[i0],fragments[i0+1],fragments[i0],fragments[l0+1],fragments[l0],n1,n2,group[i0],jsdiv);
		group1[l0+2*(y-1)]=group1[l0];gp[l0+2*(y-1)]=gp[l0];

		for (i1=0;i1<ttl_oligos;i1++)
		{
  		mat[i1][l0+2*(y-1)]=mat[i1][l0];
		}

		fragments[l0+2*(y-1)]=fragments[l0];
		fragments[l0+2*(y-1)+1]=fragments[l0+1];
		group[l0+2*(y-1)]=group[l0];
   }
}
else
   {z1++;
		//fprintf (fnonconlist,"%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%f\n",j5,group1[i0],fragments[i0+1],fragments[i0],fragments[l0+1],fragments[l0],n1,n2,group[i0],jsdiv);
		group1[l0+2*(y-1)]=group1[l0];gp[l0+2*(y-1)]=gp[l0];
  		for (i1=0;i1<ttl_oligos;i1++)
  		{
		mat[i1][l0+2*(y-1)]=mat[i1][l0];
		}

		fragments[l0+2*(y-1)]=fragments[l0];
		fragments[l0+2*(y-1)+1]=fragments[l0+1];
		group[l0+2*(y-1)]=group[l0];
		
    }
 }
		h3=h3+2*(y-1);
}
		//printf ("%d\t%d\t%d\n",h3,h,count);
   	    if(z>0)
       	     {	
		Noncontiguous (group,group1,mat,h3,gp,gp1,R,i);
       	     }

}
/*###########################################################################################################################################################*/

double entropy (double *freq_oligos,double seqLen)
{
	int j1,j2=0;
	double ent=0.0,*prob_oligos=malloc(256*sizeof(double)),sumoligos=0;	
		
	for (j1=0;j1<ttl_oligos;j1++)
	{
		sumoligos=(freq_oligos[j2]+freq_oligos[j2+1]+freq_oligos[j2+2]+freq_oligos[j2+3]);
 		if (sumoligos>0)
 		{
			prob_oligos[j1]=(double)freq_oligos[j1]/sumoligos;
			
			if (prob_oligos[j1]>0.0)
				ent+=freq_oligos[j1]*(log (prob_oligos[j1]));
 		}
		
		if (((j1+1)%4 ==0)&&(j1!=0))
			{j2+=4;/*printf ("\n");*/}
	}
	free(prob_oligos);
	return ent;
}

/**********************************************************************/
float gammp(float a, float x)
	{
		void gcf(float *gammcf, float a, float x, float *gln);
		void gser(float *gamser, float a, float x, float *gln);
		void nrerror(char error_text[]);
		float gamser, gammcf,gln;

		if(x < 0.0 || a <= 0.0) nrerror("invalid arguments in routine gammp ");
		if (x< (a+1.0)){
			gser(&gamser, a,x,&gln);
			return gamser;
			}else {
			gcf(&gammcf,a,x,&gln);
			return 1.0-gammcf;
			}
	}

/**********************************************************************/
void gser(float *gamser, float a,float x, float *gln)
	{
	//float gammln(float xx);
	void nrerror(char error_text[]);
	int n;
	float sum,del,ap;
	*gln=gammln(a);
	if (x <= 0.0){
	if (x< 0.0) nrerror ("x less than 0 in routine gser");
	*gamser=0.0;
	return;
	}else{
	ap=a;
	del=sum=1.0/a;
	for(n=1;n<=ITMAX;n++){
	++ap;
	del*=x/ap;
	sum+=del;
	if(fabs(del)<fabs(sum)*EPS) {
		*gamser=sum*exp(-x+a*log(x)-(*gln));
		return;
		}
	}
	nrerror("a too large,ITMAX too small in routine gser");
	return;
	}
	}

/**********************************************************************/
void gcf(float *gammcf, float a, float x, float *gln)
	{
		//float gammln (float xx);
		void nrerror (char error_text[]);
		int i;
		float an,b,c,d,del,h;
		*gln=gammln(a);
		b=x+1.0-a;
		c=1.0/FPMIN;
		d=1.0/b;
		h=d;
		for(i=1;i<=ITMAX;i++)
		{
			an=-i*(i-a);
			b+=2.0;
			d=an*d+b;
			if(fabs(d)<FPMIN) d=FPMIN;
			c=b+an/c;
			if(fabs(c)< FPMIN) c=FPMIN;
			d=1.0/d;
			del=d*c;
			h*=del;
			if (fabs(del-1.0)<EPS) break;
		}

			if(i> ITMAX) nrerror("a too large,ITMAX too small in gcf");
			*gammcf=exp(-x+a*log(x)-(*gln))*h;
	}

/**********************************************************************/
float gammln(float xx)
	{
		double x,y,tmp,ser;
		static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091, -1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
		int j;
		y=x=xx;
		tmp=x+5.5;
		tmp-=(x+0.5)*log(tmp);
		ser=1.000000000190015;
		for(j=0;j<=5;j++) ser+= cof[j]/++y;
		return -tmp+log(2.5066282746310005*ser/x);
	}
/**********************************************************************/
void nrerror (char error_text[])
	{
	printf("%s\n",error_text);
	}
/***************************************************************************************************************************************************************************************************/

	
