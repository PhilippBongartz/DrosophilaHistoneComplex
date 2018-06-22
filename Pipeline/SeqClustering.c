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

#define maxlength2 60000   // max read length
#define maxlength1 60000   // max template piece length


//#define PRINT

/*
This is a primitive dna sequence clustering tool.
We start by calculating the alignment scores of random pairs.
Then we build clusters by connecting sequences that share overlaps.
*/




void *Guarded_Malloc(size_t size)
{ void *p;

  p = malloc(size);
  if (p == NULL)
    { fprintf(stderr,"\nError in SeqClustering\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}

void *Guarded_Realloc(void *p, size_t size)
{ p = realloc(p,size);
  if (p == NULL)
    { fprintf(stderr,"\nError in SeqClustering:\n");
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

//ReadingFasta -> Reads: With offset, so that one read after another can be read.
off_t offset=0; //global
int readcount=-1;
char *Read;
int readlength=0;
char *buffer;
int readanzahl;
int ReadingFasta(FILE * File)
{

  char *s;
  //FILE * File=NULL;
  size_t len=70000;

  //File=fopen(Path,"r");
  if(fseek(File,offset,SEEK_SET)){printf("Problem.\n");exit(1);}

  int i;

  //printf("%p\n",File );
  fflush(stdout);

  if(File==NULL){exit(1);}

  //File+=offset;
  int bytes=0;
  int basecount=0;
  int nextread=0;
  while(nextread<2)
  {
    s = fgets(buffer, len, File);
    if(s==NULL){offset=0;return 1;}

    if(buffer[0]=='>')
    {
      nextread+=1;
      i=0;
      while(buffer[i]!='\n' && nextread<2)  //The start of the read
      {
        i+=1;
      }
      bytes+=i+1;  //+1 because of \n
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
      bytes+=i+1;  //+1 because of  \n
    }

  }
  readcount++;

  //fclose(File);
  offset+=bytes-1;  //-1 to go back to '>'
  //printf("offset: %lld\n",offset );
  Read[basecount]='\n';

  readlength=basecount;

  #ifdef PRINT
  printf("Readlength %d\n",basecount );
  printf("Read %d\n",readcount );
  fflush(stdout);
  #endif

  return 0;
}

//Reading in the second read
char *Read2;
int readlength2=0;
int ReadingFasta2(FILE * File)
{
  char *s;
  size_t len=70000;

  if(fseek(File,offset,SEEK_SET)){printf("Problem.\n");exit(1);}

  int i;

  //printf("%p\n",File );
  fflush(stdout);

  if(File==NULL){exit(1);}

  //File+=offset;
  int bytes=0;
  int basecount=0;
  int nextread=0;

  while(nextread<2)
  {
    s = fgets(buffer, len, File);
    if(s==NULL){offset=0;return 1;}

    if(buffer[0]=='>')
    {
      nextread+=1;
      i=0;
      while(buffer[i]!='\n' && nextread<2)  //Start of the read
      {
        i+=1;
      }
      bytes+=i+1;  //+1 because of \n
    }
    else if(nextread==1)  //The bases are read in 
    {
      //printf("\n\n%d\n",readcount);
      i=0;
      while(buffer[i]!='\n')
      {
        if(buffer[i]=='A' || buffer[i]=='a')Read2[basecount]='a';
        else if(buffer[i]=='C' || buffer[i]=='c')Read2[basecount]='c';
        else if(buffer[i]=='G' || buffer[i]=='g')Read2[basecount]='g';
        else if(buffer[i]=='T' || buffer[i]=='t')Read2[basecount]='t';
        else {basecount--;}        
        basecount++;
        //printf("%c ",buffer[i]);
        i++;
      }
      bytes+=i+1;  //+1 because of \n
    }

  }
  readcount++;

  //printf("offset: %lld\n",offset );
  Read2[basecount]='\n';
  readlength2=basecount;

  #ifdef PRINT
  printf("Readlength2 %d\n",basecount );
  printf("Read2 %d\n",readcount );
  fflush(stdout);
  #endif

  return 0;
}

//Inverting reads
char BaseInverse[300];
void Read2Inversion()
{
  int i;
  char temp;
  for(i=readlength2/2;i>-1;i--)
  {
    temp=Read2[i];
    Read2[i]=BaseInverse[(int)Read2[readlength2-i-1]];
    Read2[readlength2-i-1]=BaseInverse[(int)temp];
  }
}


//First we need an overlapping aligning function
int **Matrix;    //The matrix except first row and left column, those are used to initialise edge cases, accessed by -1
int **Schatten1;   //The whole matrix
int **Schatten2;   //Except left column
int maxlength;
float Aligner(int minimaloverlap)
{
  //int length1=strlen(shortstring);
  //int length2=strlen(longstring);

  int x,y;

  //printf("%d %d %d %d\n",minimaloverlap, readlength, readlength2, maxlength);fflush(stdout);
  //Initialising:
  for(x=-1;x<readlength;x++)Matrix[x][-1]=0;
  for(y=0;y<readlength2;y++)Matrix[-1][y]=0;   

  //Killing overlaps that are too short:
  for(x=readlength-1;x>readlength-minimaloverlap;x--)Matrix[x][-1]=readlength+readlength2;
  for(y=readlength2-1;y>readlength2-minimaloverlap;y--)Matrix[-1][y]=readlength+readlength2;  


  //Matrix filling
  int m;
  for(x=0;x<readlength;x++)
  {
    for(y=0;y<readlength2;y++)
    {
      m=1;  
      if(Read[x]==Read2[y])m=0;
      Matrix[x][y]=Matrix[x-1][y-1]+m;
      Matrix[x][y]=intmin(Matrix[x-1][y]+1,Matrix[x][y]);
      Matrix[x][y]=intmin(Matrix[x][y-1]+1,Matrix[x][y]);
    }
  }

  //Backtracking
  x=readlength-1;
  y=readlength2-1;

  float min;
  int i,einstieg_y,einstieg_x;
  einstieg_y=y;
  einstieg_x=x;
  min=(float)Matrix[einstieg_x][einstieg_y]/(float)(intmin(readlength,einstieg_y));

  for(i=readlength2-1;i>1;i--)
  { 
    if((float)Matrix[x][i]/(float)i <min)
    {
      min=(float)Matrix[x][i]/(float)(intmin(readlength,i));
      einstieg_y=i;
    }
  }
  
  for(i=readlength-1;i>1;i--)
  { 
    if((float)Matrix[i][y]/(float)(intmin(readlength2,i))<min)
    {
      min=(float)Matrix[i][y]/(float)(intmin(readlength2,i));
      einstieg_x=i;
      einstieg_y=y;  //superseded
    }
  }
  //printf("Overlap length %d\n",intmin(einstieg_x, einstieg_y) );
  //printf("%d/%d %d/%d\n",einstieg_x,readlength,einstieg_y,readlength2 );fflush(stdout);
  return min;
}

int *offsets;
int *SeqClassArray;
int ReadCounter(FILE * File, int nonrepnumber)
{
  char *s;
  //FILE * File=NULL;
  size_t len=70000;
  //File=fopen(Path,"r");
  int count=0;
  if(File==NULL){exit(1);}
  int i;
  int length=0;
  int bytes=0;
  int nonrepcount=0;
  maxlength=0;

  offsets=Guarded_Malloc(sizeof(int)*nonrepnumber);
  while((s = fgets(buffer, len, File)) != NULL)
  {
    if(buffer[0]=='>')
    {
      if(SeqClassArray[count]) //offsets of the non-repetitive sequences
      {
        //if(length>maxlength){maxlength=length;}
        offsets[nonrepcount]=bytes-1;
        nonrepcount++;
      }

      if(count>0 && SeqClassArray[count-1])  //the length is of the predecessor read
      {
        if(length>maxlength){maxlength=length;}
      }    

      count++;
      length=0;
      i=0;
      while(buffer[i]!='\n'){bytes++;i++;}
    }
    else
    {
      i=0;
      while(buffer[i]!='\n')
      {
        i++;
        length++;
        bytes++;
      }
    }
    bytes++; // \n
  }
  if(length>maxlength){maxlength=length;}
  //fclose(File);

  return count;
}

//Reading in the Seqclass information: Which seqs are non-repetitive. 
int SeqClass(char *Path)
{
  char *s;
  FILE * File=NULL;
  size_t len=70000;
  File=fopen(Path,"r");
  int count=0;
  int yes=0;
  if(File==NULL){exit(1);}
  int estimate=20000;
  buffer=Guarded_Malloc(sizeof(char)*70000);
  SeqClassArray=Guarded_Malloc(sizeof(int)*estimate);

  while((s = fgets(buffer, len, File)) != NULL)
  {
    if(buffer[0]=='l'){SeqClassArray[count]=1;yes++;}
    if(buffer[0]=='r')SeqClassArray[count]=0;
    count++;
    if(count>estimate-1)
    {
      estimate+=1000;
      SeqClassArray=Guarded_Realloc(SeqClassArray,sizeof(int)*estimate);
    }
  }
  fclose(File);
  //printf("read number according to seqclass %d\n",count );
  return yes;
}

//Memory allocation
void MatrixMemory(int maxlength)
{
  Schatten1=Guarded_Malloc(sizeof(int*)*(maxlength+2));
  int x;
  for(x=0;x<maxlength+2;x++)*(Schatten1+x)=Guarded_Malloc(sizeof(int)*(maxlength+2));
  Schatten2=Guarded_Malloc(sizeof(int*)*(maxlength+2));
  for(x=0;x<maxlength+2;x++)*(Schatten2+x)=(*(Schatten1+x)+1); 
  Matrix=Schatten2+1;  

  BaseInverse['a']='t';
  BaseInverse['t']='a';
  BaseInverse['c']='g';
  BaseInverse['g']='c';
  BaseInverse['A']='T';
  BaseInverse['T']='A';
  BaseInverse['C']='G';
  BaseInverse['G']='C';

  Read=Guarded_Malloc(sizeof(char)*(maxlength+1));
  Read2=Guarded_Malloc(sizeof(char)*(maxlength+1));

}

//This function groups n-connected seqs into clusters.
int* SeqGrouping(double** ScoreMatrix, int n, double error_cutoff, int counti)
{
  int *doubleconnection=Guarded_Malloc(sizeof(int)*counti);
  int *cluster_no=Guarded_Malloc(sizeof(int)*counti);
  int i,j,l;
  for(i=0;i<counti;i++)cluster_no[i]=i;
  int mini_no;
  int temp_no;

  //For each seq we look at the connections.
  for(i=0;i<counti;i++)
  {
    for(j=0;j<counti;j++)doubleconnection[j]=0;
    for(j=0;j<counti;j++)
    {
      //For each connection we collect the connections.
      if(ScoreMatrix[i][j]<error_cutoff)
      {
        for(l=0;l<counti;l++)
        {
          if(ScoreMatrix[l][j]<error_cutoff)
          {
            //We count which of the double connections occur n times. 
            doubleconnection[l]++;  //double connection from i to l
          }
        }
      }
    }
    //The minimum cluster_no:
    mini_no=cluster_no[i];
    for(l=0;l<counti;l++)
    {
      if(doubleconnection[l]>=n)
      {
        if(cluster_no[l]<mini_no)
        {
          mini_no=cluster_no[l];
        }
      }
    }
    //All seqs with the cluster_nos involved are changed to mini_no
    for(l=0;l<counti;l++)
    {
      if(doubleconnection[l]>=n || l==i)
      {
        temp_no=cluster_no[l];
        for(j=0;j<counti;j++)
        {
          if(cluster_no[j]==temp_no)
          {
            cluster_no[j]=mini_no;
          }
        }
      }
    }
  }

  printf("Cluster numbers: ");
  for(i=0;i<counti;i++)printf("%d ",cluster_no[i]);
  printf("\n");

  return cluster_no;
}

// outputting the clustering result.
void Output(char *outputfile, int* cluster_no, int counti)
{
  // The indices of the unique sequences:
  int *seq_no=Guarded_Malloc(sizeof(int)*counti);
  int i;
  int j=0;
  for(i=0;i<readanzahl;i++)
  {
    if(SeqClassArray[i])
    {
      seq_no[j]=i;
      j++;
    }
  }
  // The cluster_no with more than one member:
  int *no_count=Guarded_Malloc(sizeof(int)*counti);
  for(i=0;i<counti;i++)no_count[i]=0;
  for(i=0;i<counti;i++)no_count[cluster_no[i]]++;
  // Outputting those clusters unique sequence indices:
  FILE *datei;
  datei=fopen(outputfile,"w");
  for(i=0;i<counti;i++)
  {
    if(no_count[i]>1)
    {
      for(j=0;j<counti;j++)
      {
        if(cluster_no[j]==i)
        {
          fprintf(datei,"%d ",seq_no[j]);
        }
      }
      fprintf(datei,"\n");
    }
  }
  fclose(datei);
}

void Help()
{
  printf("Usage: ./SeqClustering Seq.fasta SeqClass\n");
  printf("Flags:\n");
  printf("-o clustering_path    Path of the resulting clustering information. Default: SimulatedSeqClustering.\n");
  printf("-e <0.30>             The alignment error cutoff being used to sort seqs into the same cluster.\n");
  printf("-m <500>              The minimal overlap length.\n");
  printf("-n <2>                The minimal number of shared overlaps to be assigned to the same cluster.\n");
  exit(0);
}




int main(int argc, char *argv[])
{

  char *Readspath_p;
  char *SeqClasspath_p;
  if(argc<2){printf("Usage: ./SeqClustering Seq.fasta SeqClass\n");exit(0);}
  Readspath_p=argv[1];
  SeqClasspath_p=argv[2];

  FILE * File=NULL;
  File=fopen(Readspath_p,"r");
  if(File==NULL){printf("%s could not be loaded.\n",SeqClasspath_p); exit(1);}

  //Specifying the output path: 
  char outputfile[]="SeqClusters";
  char *output_p=&outputfile[0];

  int i,j;

  //arguments: -o outputfile, -h help, , -e error_cutoff,

  double error_cutoff=0.30;
  int minimaloverlap=500;
  int n=2;

  for(i=1;i<argc;i++)
  {

    if(argv[i][0]=='-' && argv[i][1]=='o')
    {
      //printf("%s\n",argv[i]);
      if(i+1<argc)output_p=argv[i+1];
    }     

    if(argv[i][0]=='-' && argv[i][1]=='e')
    {
      error_cutoff=atof(argv[i+1]);
    }   

    if(argv[i][0]=='-' && argv[i][1]=='m')
    {
      minimaloverlap=atoi(argv[i+1]);
    } 

    if(argv[i][0]=='-' && argv[i][1]=='n')
    {
      n=atoi(argv[i+1]);
    } 

    if(argv[i][0]=='-' && argv[i][1]=='h')
    {
      Help();
    }    
  }

  // First we determine which sequences are non-repetitive - those will be clustered:
  int counti=SeqClass(SeqClasspath_p);

  //Initialising the ScoreMatrix
  double **ScoreMatrix=Guarded_Malloc(sizeof(double*)*counti);
  for(i=0;i<counti;i++)ScoreMatrix[i]=Guarded_Malloc(sizeof(double)*counti);
  for(i=0;i<counti;i++)
  {
    for(j=0;j<counti;j++)
    {  
      ScoreMatrix[i][j]=1.0;
    }
  }
  for(i=0;i<counti;i++)ScoreMatrix[i][i]=0.0;

  // Then we find out how to access exactly these non-repetitive sequences:
  readanzahl=ReadCounter(File,counti);

  printf("Read count %d, maxlength %d\n",readanzahl,maxlength);fflush(stdout);
  printf("Number of sequences: %d\n",counti );


  //Memory allocation:
  MatrixMemory(maxlength);
  printf("Memory allocated.\n");

  //ReadingFasta(Readspath_p);  
  int histo[100];
  for(i=0;i<100;i++)histo[i]=0;

  int prozent=5;
  for(i=0;i<counti;i++)
  {
    offset=offsets[i];
    ReadingFasta(File);
    for(j=i+1;j<counti;j++)
    {
      if(readlength>minimaloverlap)
      {
        offset=offsets[j];
        ReadingFasta2(File);

        if(readlength2>minimaloverlap)
        {
          float score=Aligner(minimaloverlap);
          ScoreMatrix[i][j]=score;
          ScoreMatrix[j][i]=score;
          //printf("Score %f\n",score );
          if(score>error_cutoff) //Inverting the read and aligning again
          {
            Read2Inversion();
            float score=Aligner(minimaloverlap);
            if(score<ScoreMatrix[i][j])
            {
              ScoreMatrix[i][j]=score;
              ScoreMatrix[j][i]=score;            
            }
          }
          histo[(int)(ScoreMatrix[i][j]*100)]++;
        }
      }
    }

    //SeqGrouping(ScoreMatrix, n, error_cutoff, counti);
    //printf("\n");
    //printf("Alignment score histogramm:\n");
    //for(j=0;j<100;j++)printf("%d ",histo[j]);
    //printf("\n\n");

    //if(i*100/readanzahl>prozent)  //

    if((counti*i-(i*i)/2)*100/((counti*counti)/2)>prozent)
    {
      printf("%d %% done.\n",prozent );
      prozent+=5;
    }
  }

  int *cluster_no=SeqGrouping(ScoreMatrix, n, error_cutoff, counti);

  fclose(File);
  //Output: 
  printf("Outputting results.\n");
  Output(output_p, cluster_no, counti);


  exit(0);


} 





