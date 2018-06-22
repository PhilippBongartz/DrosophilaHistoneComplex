#include <stdio.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>

#include <stdint.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/stat.h>

#define maxlength2 50000
#define maxlength1 50000

//#define PRINT



/*
This tool creates an initial multiple sequence alignment by aligning reads to a template.
The resulting msa is only a rough arrangement of the most important bases and has to be refined.
*/


void *Guarded_Malloc(size_t size)
{ void *p;

  p = malloc(size);
  if (p == NULL)
    { fprintf(stderr,"\nError in InitialAligner\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}

void *Guarded_Realloc(void *p, size_t size)
{ p = realloc(p,size);
  if (p == NULL)
    { fprintf(stderr,"\nError in InitialAligner:\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}

int intmin(int x, int y){if(x<y){return x;} return y;}

int intmax(int x, int y){if(x>y){return x;} return y;}

int intabs(int a)
{
  if(a<0){return -a;}
  return a;
}


//ReadingFasta -> Reads: Reads in read by read.
off_t offset=0; //global
int readcount=0;
char Read[70000];
int readlength=0;
int readanzahl;
char buffer[70000];
int ReadingFasta(char *Path)
{
  char *s;
  FILE * File=NULL;
  size_t len=70000;

  File=fopen(Path,"r");
  if(fseek(File,offset,SEEK_SET)){printf("The End.\n");fflush(stdout);exit(1);}

  int i;

  //printf("%p\n",File );
  fflush(stdout);

  if(File==NULL){printf("File == NULL.\n");fflush(stdout);exit(1);}

  //File+=offset;
  int bytes=0;
  int basecount=0;
  int nextread=0;
  while(nextread<2)
  {
    s = fgets(buffer, len, File);
    if(s==NULL)
    {
    	if(feof(File))nextread+=1;
    	else{return 1;}
    }

    else if(buffer[0]=='>')
    {
      nextread+=1;
      i=0;
      while(buffer[i]!='\n' && nextread<2)  //The beginning of the read
      {
        i+=1;
      }
      bytes+=i+1;  //+1 wegen \n
    }
    else  //The bases are read in 
    {
      //printf("\n\n%d\n",readcount);
      i=0;
      while(buffer[i]!='\n')
      {
        if(buffer[i]=='A' || buffer[i]=='a')Read[basecount]='a';
        else if(buffer[i]=='C' || buffer[i]=='c')Read[basecount]='c';
        else if(buffer[i]=='G' || buffer[i]=='g')Read[basecount]='g';
        else if(buffer[i]=='T' || buffer[i]=='t')Read[basecount]='t';
        else {basecount--;}        
        basecount++;
        //printf("%c ",buffer[i]);
        i++;
      }
      bytes+=i+1;  //+1 wegen \n
    }

  }

  fclose(File);
  offset+=bytes-1;  //-1 to go back to '>'
  //printf("offset: %lld\n",offset );
  Read[basecount]='\n';
  readlength=basecount;

  #ifdef PRINT
  printf("Read %d, Readlength %d\n",readcount,basecount );
  #endif

  readcount++;
  return 0;
}


//Reading in the template
char Template[70000];
int templatelength=0;
void ReadingTemplate(char *Path)
{
  char *s;
  FILE * File=NULL;
  size_t len=70000;
  File=fopen(Path,"r");

  int i;

  if(File==NULL){exit(1);}

  int basecount=0;

  while((s = fgets(buffer, len, File)) != NULL)
  {
    if(buffer[0]!='>')
    {
      i=0;
      while(buffer[i]!='\n')
      {
        if(buffer[i]=='A' || buffer[i]=='a')Template[basecount]='a';
        else if(buffer[i]=='C' || buffer[i]=='c')Template[basecount]='c';
        else if(buffer[i]=='G' || buffer[i]=='g')Template[basecount]='g';
        else if(buffer[i]=='T' || buffer[i]=='t')Template[basecount]='t';
        else {basecount--;}
        basecount++;
        i++;
      }
    }

  }

  fclose(File);
  Template[basecount]='\n';
  templatelength=basecount;

}


int **Matrix;    //The matrix except left column and upper row, those are accessed via -1 for edge case initialisation.
int **Schatten1;   //The full matrix
int **Schatten2;   //The matrix except the left column
int **Alignments;
int *Readlengths;
double *AlignmentError;
int aligncount;
int IntoAligner(char *shortstring, char *longstring, int length1, int length2)
{
  //int length1=strlen(shortstring);
  //int length2=strlen(longstring);

  int x,y;


  //Initialisation
  for(x=-1;x<length1;x++)Matrix[x][-1]=x+1;
  for(y=0;y<length2;y++)Matrix[-1][y]=0;    //Seq1 is aligned into Seq2

  //Matrixfüllen
  int m;
  for(x=0;x<length1;x++)
  {
  	for(y=0;y<length2;y++)
  	{
  	  m=1;	
  	  if(shortstring[x]==longstring[y])m=0;
  	  Matrix[x][y]=Matrix[x-1][y-1]+m;
      Matrix[x][y]=intmin(Matrix[x-1][y]+1,Matrix[x][y]);
      Matrix[x][y]=intmin(Matrix[x][y-1]+1,Matrix[x][y]);
  	}
  }

  //Backtracking
  char *EditScript;
  int count=0;
  EditScript=Guarded_Malloc(sizeof(char)*(length1+length2));
  x=length1-1;
  y=length2-1;

  int min,i,einstieg_y;
  einstieg_y=y;
  min=Matrix[x][einstieg_y];
  for(i=length2-1;i>0;i--)
  {	
  	if(Matrix[x][i]<min)
  	{
  		min=Matrix[x][i];
  		einstieg_y=i;
  	}
  }

  #ifdef PRINT
  printf("Entry %d\n",einstieg_y );
  printf("%d\n",Matrix[length1-1][einstieg_y]);
  #endif

  AlignmentError[aligncount]=(double)Matrix[length1-1][einstieg_y]/(double)length1;

  #ifdef PRINT
  printf("Alignment Error[%d]=%f\n",aligncount,AlignmentError[aligncount] );
  #endif

  while(y>einstieg_y)
  {
    EditScript[count]='i';  
    count++;  
    y--;  	
  }


  while(x>-1 && y>-1)
  {
  	//printf("%d %d\n",x,y );
  	//printf("%d %d %d -> %d\n",Matrix[x-1][y-1],Matrix[x][y-1],Matrix[x-1][y],Matrix[x][y] );
  	m=1;	
  	if(shortstring[x]==longstring[y])m=0;

  	if(Matrix[x][y]==Matrix[x-1][y-1]+m)
  	{
  	  if(m) { EditScript[count]='s';}
  	  else  {EditScript[count]='m';}
  	  count++;
  	  x--;
  	  y--;
  	}
    else if(Matrix[x][y]==Matrix[x][y-1]+1)
    {
      EditScript[count]='i';  
      count++;  
      y--;
    }    
  	else if(Matrix[x][y]==Matrix[x-1][y]+1)
  	{
  	  EditScript[count]='d';	
  	  count++;
  	  x--;  		
  	}
  	else
  	{
  	  printf("Backtracking error.\n");
  	  exit(1);
  	}	
  	//printf("x:%d y:%d Matr:%d\n",x,y,Matrix[x][y]);
  }

  while(x>-1)
  {
  	EditScript[count]='d';	
    count++;
		x--;  	  	
  } 
  while(y>-1) 	
  {
  	EditScript[count]='i';	
  	count++;	
  	y--;  	
  }

  //Inverting edit script:
  for(x=0;x<count/2;x++)
  {
  	m=EditScript[x];
  	EditScript[x]=EditScript[count-x-1];
  	EditScript[count-x-1]=m;
  }

/*  for(x=0;x<count;x++)printf("%c",EditScript[x]);
  printf("\n");	*/

  // Turning the edit script into a alignment:

	Alignments[aligncount]=Guarded_Malloc(sizeof(int)*length1);
	Readlengths[aligncount]=length1;

  x=0;
  y=0;
	for(i=0;i<count;i++)
	{
		if(EditScript[i]=='s' || EditScript[i]=='m')
		{
			//printf("%c - %c\n", shortstring[x],longstring[y]);
			Alignments[aligncount][x]=y;
			x++;
			y++;
		}
		if(EditScript[i]=='i')
		{
			//printf("_ - %c\n", longstring[y]);		
			y++;	
		}
		if(EditScript[i]=='d')
		{
			//printf("%c - _\n", shortstring[x]);			
			Alignments[aligncount][x]=-1;
			x++;
		}
	}
/*	for(x=0;x<length1;x++)printf("%d ",Alignments[aligncount][x] );
	printf("\n");*/

	aligncount++;


  return min;
}

int ReadCounter(char *Path)
{
  char *s;
  FILE * File=NULL;
  size_t len=70000;
  File=fopen(Path,"r");
  int count=0;
  if(File==NULL){exit(1);}
  while((s = fgets(buffer, len, File)) != NULL)
  {
    if(buffer[0]=='>')count++;
	}
	fclose(File);
	return count;
}

//Memory allocation
void MatrixMemory()
{
  Schatten1=Guarded_Malloc(sizeof(int*)*(maxlength1+1));
  int x;
  for(x=0;x<maxlength1+1;x++)*(Schatten1+x)=Guarded_Malloc(sizeof(int)*(maxlength2+1));
  Schatten2=Guarded_Malloc(sizeof(int*)*(maxlength1+1));
  for(x=0;x<maxlength1+1;x++)*(Schatten2+x)=(*(Schatten1+x)+1); 
  Matrix=Schatten2+1;  

	Alignments=Guarded_Malloc(sizeof(int*)*readanzahl);
	Readlengths=Guarded_Malloc(sizeof(int)*readanzahl);
  AlignmentError=Guarded_Malloc(sizeof(double)*readanzahl);
}


void Building_MSA(char *outputfile, char *seqclassfile, char *Readspath_p, double errorcutoff)
{
  FILE *datei;
  datei=fopen(outputfile,"w");


  FILE *datei2;
  datei2=fopen(seqclassfile,"w");
  //Alignment[j][i] tells us where the i-th base of the j-th read went. 
  // ==-1 indicates that the i-th base is placed between two template bases.

	//Determine the maximal gap between all bases:
	int *Gapcount=Guarded_Malloc(sizeof(int)*templatelength+1);
	int i,j,k,l,count,gap;
	for(i=0;i<templatelength+1;i++)Gapcount[i]=0;

	for(j=0;j<readanzahl;j++)
	{
		gap=0;
		count=0;
		//first gap 
		i=0;
		while(Alignments[j][i]==-1)i++;
		gap=Alignments[j][i];  //The bases before 'gap'

		for(i=0;i<Readlengths[j];i++)
		{
			if(Alignments[j][i]==-1)
			{
				count++;
				if(count>Gapcount[gap]){Gapcount[gap]=count;}
			}
			else
			{
				gap=Alignments[j][i]+1;  //The bases before Alignments[j][i] are already counted.
				count=0;
			}
		}
	}

	//Write out the reads while respecting the gaps:
	for(j=0;j<readanzahl;j++)
	{
		ReadingFasta(Readspath_p);
    if(AlignmentError[j]<errorcutoff)
    {
      fprintf(datei2,"r\n");
  		if(readlength>0)
  		{
  			k=0;
  			i=0;
  			while(i<templatelength+1)
  			{
  				//gap bases into the gap
  				count=0;
  				while(k<readlength && Alignments[j][k]==-1)
  				{
  					fprintf(datei,"%c",Read[k]);
  					if(Read[k]=='\n'){printf("Bumm\n"); exit(0);}
  					k++;
  					count++;
  				}
  				//Filling the gap with '-'
  				for(l=count;l<Gapcount[i];l++)fprintf(datei,"-");

  				//Next base base if it exists
  				if(k<readlength && Alignments[j][k]==i)
  				{
  					fprintf(datei,"%c",Read[k]);
  					if(Read[k]=='\n'){printf("Bumm\n"); exit(0);}
  					k++;
  				}
  				else
  				{
  					fprintf(datei,"-");
  				}
  				i++;
  			}
  		}
  		else
  		{
  			for(i=0;i<templatelength+1;i++)
  			{
  				fprintf(datei,"-");
  				for(k=0;k<Gapcount[i];k++)fprintf(datei,"-");
  			}
  		}
  		fprintf(datei,"\n"); 
    }
		//if(j==410)for(i=0;i<readlength;i++)printf("%d\n", Alignments[j][i]);

    else //No repeat seq, 'l' is for large scale variation as opposed to 'r' for repetitive. 
    {
      fprintf(datei2,"l\n");
    }
	}
	fclose(datei);
  fclose(datei2);
}


void Help()
{
  printf("Usage: ./InitialAligner template.fasta Seq.fasta\n");
  printf("Flags:\n");
  printf("-o msa_path    Path of the resulting multiple sequence alignment. Default: SimulatedMSA.\n");
  printf("-s <150>       Path of the seq class information: Which seqs are instances of the repeat. Default: SimulatedSeqClass\n");
  printf("-e <0.30>      The mapping error cutoff being used to detect instances of the template.\n");
  exit(0);
}









int main(int argc, char *argv[])
{

  char *Templatepath_p;
  char *Readspath_p;
  if(argc<2){printf("Usage: ./InitialAligner template.fasta Seq.fasta\n");exit(0);}
  Readspath_p=argv[2];
  Templatepath_p=argv[1];


  ReadingTemplate(Templatepath_p);
  printf("template length %d\n",templatelength );

  //Hier kann man das outputfile angeben: 
  char outputfile[]="SimulatedMSA";
  char *output_p=&outputfile[0];

  char seqclassfile[]="SimulatedSeqClass";
  char *seqclass_p=&seqclassfile[0];

  double errorcutoff=0.30;

  int i;

  for(i=1;i<argc;i++)
  {

    if(argv[i][0]=='-' && argv[i][1]=='o')
    {
      //printf("%s\n",argv[i]);
      if(i+1<argc)output_p=argv[i+1];
    }  

    if(argv[i][0]=='-' && argv[i][1]=='s')
    {
      //printf("%s\n",argv[i]);
      if(i+1<argc)seqclass_p=argv[i+1];
    }  

    if(argv[i][0]=='-' && argv[i][1]=='e')
    {
      errorcutoff=atof(argv[i+1]);
    }   

    if(argv[i][0]=='-' && argv[i][1]=='h')
    {
      Help();
    } 
  }

  printf("output file: %s\n",output_p );

  readanzahl=ReadCounter(Readspath_p);
  printf("read count %d\n",readanzahl );

  //readanzahl=100;

  //Memory allocation:
  MatrixMemory();
  aligncount=0;
  readcount=0;
  int prozent=5;

  for(i=0;i<readanzahl;i++)
  {
	  ReadingFasta(Readspath_p);
	  IntoAligner(Read, Template, readlength, templatelength);

    if(i*100/readanzahl>prozent)
    {
      printf("%d %% done.\n",prozent );
      prozent+=5;
    }
  }

	offset=0; //global
	readcount=0;
  printf("Writing the msa.\n");
	Building_MSA(output_p, seqclass_p, Readspath_p, errorcutoff);

  exit(0);
}  

