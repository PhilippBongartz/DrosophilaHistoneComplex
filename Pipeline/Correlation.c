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

//#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>

#define Max_Var_Anzahl 200000
#define Max_Sig_Anzahl 20000


/*
This tool calculates the statistical significance of intersections
between basegroups of different columns of a multiple sequence alignment.
It outputs signatures build from the bases of columns with strong correlations.
-b calculates only the correlations of columns with a majority of bases
-c x uses the correlation logprob cutoff x (default 12)
-o file specifies the output file (default "Signatures")
-n x uses the x best columns above the cutoff to build the signatures
-v verbose

It requires gnu scientific library:
gcc -Wall -I/usr/local/include -c Correlation.c 
gcc -L/usr/local/lib Correlation.o -lgsl -lgslcblas -lm -o Correlation
*/



unsigned long *Groups[Max_Var_Anzahl*5];
unsigned long *LocalCoverage[Max_Var_Anzahl];
int v;

void *Guarded_Malloc(size_t size)
{ void *p;

  p = malloc(size);
  if (p == NULL)
    { fprintf(stderr,"\nError in Correlation\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}

void *Guarded_Realloc(void *p, size_t size)
{ p = realloc(p,size);
  if (p == NULL)
    { fprintf(stderr,"\nError in Correlation:\n");
      fprintf(stderr,"   Out of memory\n");
      exit (1);
    }
  return (p);
}


char **Signatures;

/*******************************Group-Handling**************************************/


// BitCounter taken from Hacker's Delight, Warren, Jr. Henry S.: Die const uint kommen in .h
const uint64_t m1  = 0x5555555555555555; //binary: 0101...
const uint64_t m2  = 0x3333333333333333; //binary: 00110011..
const uint64_t m4  = 0x0f0f0f0f0f0f0f0f; //binary:  4 zeros,  4 ones ...
const uint64_t m8  = 0x00ff00ff00ff00ff; //binary:  8 zeros,  8 ones ...
const uint64_t m16 = 0x0000ffff0000ffff; //binary: 16 zeros, 16 ones ...
const uint64_t m32 = 0x00000000ffffffff; //binary: 32 zeros, 32 ones
const uint64_t hff = 0xffffffffffffffff; //binary: all ones
const uint64_t h01 = 0x0101010101010101; //the sum of 256 to the power of 0,1,2,3...
 
int popcount_3(uint64_t x) 
{
  x -= (x >> 1) & m1;             //put count of each 2 bits into those 2 bits
  x = (x & m2) + ((x >> 2) & m2); //put count of each 4 bits into those 4 bits 
  x = (x + (x >> 4)) & m4;        //put count of each 8 bits into those 8 bits 
  return (x * h01)>>56;  //returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ... 
}

int intmax(int x, int y)
{
  if(x>y)return x;
  return y;
}

int sc; //This variable is (SigAnzahl/64 +1), the number of ints necessary to encode a base group:

// Intersection of two base groups
unsigned int Schnitt(unsigned long *group1, unsigned long *group2)
{  
  int zz;
  unsigned int zahl=0;
  unsigned long dings;
  for (zz=0; zz<sc; zz++)
  {
    dings=*(group1+zz) & *(group2+zz);
    zahl+=popcount_3(dings);                   
  }
  return zahl;  
}

void GrAdd(unsigned long *group, int element)
{
  unsigned long mask=1;
  int VecInd=element/64;
  mask=mask<<(element%64);
  *(group+VecInd)=*(group+VecInd) | mask;
}

int GrElement(unsigned long *group, int element)
{
  unsigned long mask=1;
  int VecInd=element/64;
  mask=mask<<(element%64);
  if(*(group+VecInd) & mask)return 1;
  return 0;
}

void GrNull(unsigned long* Group)
{
  int i;
  for(i=0;i<sc;i++)
  {
    *(Group+i)=0;
  }  
}

int Groupsize(unsigned long *group)
{
  unsigned long *castedgroup;
  int zz;
  int zahl=0;
  castedgroup=group;
  for (zz=0; zz<sc; zz++)
    {
      zahl+=popcount_3(castedgroup[zz]);               
    }
  return zahl;  
}


/******************************************Signatures-ReadIn*****************************************/
int siglength;
int signumber;
int Coverage[Max_Var_Anzahl];

void Einlesen(char *MApath)   
{ 

  char *s;  
  int i,j,k;
  char buffer[Max_Var_Anzahl];

 
  Signatures=Guarded_Malloc(sizeof(char*)*Max_Sig_Anzahl);
  for(i=0;i<Max_Sig_Anzahl;i++)
  {
    *(Signatures+i)=Guarded_Malloc(sizeof(char)*Max_Var_Anzahl);
  }

  signumber=-1;

  //while ((s=fgets(buffer, Max_Var_Anzahl-2, stdin)) != NULL)  //for piping the msa in

  FILE * File;
  size_t len=Max_Var_Anzahl-2;
  File=fopen(MApath,"r");
  if(File==NULL){printf("MA is missing.\n"); exit(1);}

  while ((s = fgets(buffer, len, File)) != NULL)  
  {
    signumber++;
    if(v)printf("%d\n",signumber);

    siglength = strlen(buffer);
    if (buffer[siglength-1] != '\n'){printf("Wrong format #%c# \n",buffer[siglength-1]);exit(1);}
    siglength--;  //length was including '\n'

    
    for(i=0;i<siglength;i++)
    {
      if(buffer[i]=='a' || buffer[i]=='A')
      {
        Signatures[signumber][i]=0;
      }
      else if(buffer[i]=='c' || buffer[i]=='C')
      {
        Signatures[signumber][i]=1;
      }
      else if(buffer[i]=='g' || buffer[i]=='G')
      {
        Signatures[signumber][i]=2;
      }
      else if(buffer[i]=='t' || buffer[i]=='T')
      {
        Signatures[signumber][i]=3;
      }
      else if(buffer[i]=='-' || buffer[i]=='_')
      {
        Signatures[signumber][i]=4;
      }
      else
      {
        Signatures[signumber][i]=5;
      }
    }
  }
  signumber++;  //Not the last index, the number
  sc=(signumber/64)+1;

  //Groups are organized into an array
  for(i=0;i<siglength*5;i++)
  {
    Groups[i]=Guarded_Malloc(sizeof(unsigned long)*sc);
    GrNull(Groups[i]);
  }  

  //Coverage for each group -> to calculate shared coverage of two groups
  for(i=0;i<siglength;i++)
  {
    LocalCoverage[i]=Guarded_Malloc(sizeof(unsigned long)*sc);
    GrNull(LocalCoverage[i]);
  }  

  for(i=0;i<siglength;i++)Coverage[i]=0;

  for(i=0;i<siglength;i++)
  {
    for(j=0;j<signumber;j++)
    {
      for(k=0;k<5;k++)
      {
        if(Signatures[j][i]==k) 
        {
          GrAdd(Groups[i*5+k],j);
        }
      }

      if(Signatures[j][i]<5)
      {
        GrAdd(LocalCoverage[i],j);
        Coverage[i]++;
      }
    }
  }

  if(v)printf("Siglength %d, signumber %d, Groups %d, sc %d ...\n",siglength,signumber,siglength*5,sc );
  fflush(stdout);
} 

int intmin(int a, int b)
{
  if(a<b)return a;
  return b;
}

// The cumulative hypergeometric probability
double CumHypGeo_Log(unsigned int schnitt, unsigned int gr1, unsigned int gr2, unsigned int cov)
{ 
  double posP=gsl_cdf_hypergeometric_P (schnitt, gr2, cov-gr2, gr1);  
  double posQ=gsl_cdf_hypergeometric_Q (schnitt-1, gr2, cov-gr2, gr1); 

  if(posP<posQ || schnitt==0)
  {
    posP=-1.0*log10(posP);
    if(isinf(posP) || posP>99)return 99.0;
    return posP;
  }
  posQ=-1.0*log10(posQ);
  if(isinf(posQ) || posQ>99)return 99.0;
  return posQ;
}

double fabs(double x){if(x>0)return x;return -1.0*x;}

void PoissonKorrs(int c, int b, char* outputfile, int cho)
{
  unsigned int i,j,cov,schnitt,gr1,gr2;
  double Z,mean,sigma;

  double *MaxCorr;
  int *Groupsizearray;

  int *Z_Histo=Guarded_Malloc(sizeof(int)*100);
  for(i=0;i<100;i++)Z_Histo[i]=0;


  int **ZZ_Histo=Guarded_Malloc(sizeof(int*)*100);
  for(i=0;i<100;i++)ZZ_Histo[i]=Guarded_Malloc(sizeof(int)*100);
  for(i=0;i<100;i++){for(j=0;j<100;j++)ZZ_Histo[i][j]=0;}  

  MaxCorr=Guarded_Malloc(sizeof(double)*siglength*5);
  for(i=0;i<siglength*5;i++)MaxCorr[i]=0;

  Groupsizearray=Guarded_Malloc(sizeof(int)*siglength*5);
  int maxcov=0;
  for(i=0;i<siglength*5;i++)
  {
    Groupsizearray[i]=Groupsize(Groups[i]);  
    if(Groupsizearray[i]>maxcov)maxcov=Groupsizearray[i];
    //printf("%d ",Groupsizearray[i]);
  }

  if(v)printf("Maxcov: %d\n",maxcov);

  int comparisons=0;

  int prozent=5;
  long unsigned siglsigl=(siglength/2)*siglength;  //max. Number of comparisons
  long unsigned progress;

  for(i=0;i<siglength*5;i++)
  {
    //printf("%d %d\n",Groupsizearray[i],Groupsizearray[j]);
    if(v && i%1000==0)printf("%d ",i );

    if(Groupsizearray[i]>100 && Coverage[i/5]>(maxcov*2)/3 && i%5!=4 && (Coverage[i/5]-Groupsizearray[(i/5)*5+4]>Coverage[i/5]/2 || 0==b)) // && Groupsizearray[i]<Coverage[i/5]/2)
    {
      for(j=i+300;j<siglength*5;j++)  
      {
        if(Groupsizearray[j]>100 && j%5!=4 && (Coverage[j/5]-Groupsizearray[(j/5)*5+4]>Coverage[j/5]/2 || 0==b))// && Groupsizearray[j]<Coverage[j/5]/2)
        {
          comparisons++;
          schnitt=Schnitt(Groups[i], Groups[j]);  

          cov=Schnitt(LocalCoverage[i/5],LocalCoverage[j/5]);
          if(cov>100)
          {
            mean=Groupsizearray[i]*Groupsizearray[j]/cov;
            sigma=sqrt(mean);
            if(fabs(schnitt-mean)/sigma>7)  //Pre-selection
            {
              gr1=Schnitt(Groups[i],LocalCoverage[j/5]);
              gr2=Schnitt(Groups[j],LocalCoverage[i/5]);
              if(gr1>0 && gr2>0)  // && (gr1<cov/2 || gr2<cov/2))
              {
                Z=CumHypGeo_Log(schnitt,gr1,gr2,cov);
                Z_Histo[(int)Z]++;
                if(Z>MaxCorr[i])MaxCorr[i]=Z;
                if(Z>MaxCorr[j])MaxCorr[j]=Z;
              }
            }
          }
        }
      }  
    }

    progress=(unsigned long)(siglength);
    progress-=(unsigned long)((i/5)/2);
    progress*=(unsigned long)(i/5);
    progress*=100;  //In percent
    //printf("%lu \n",progress/siglsigl);
    if(progress/siglsigl>prozent-1) //(siglength*5*i-(i*i)/2)*8
    {
      printf("%d %% done.\n",prozent );
      prozent+=5;
    }

  }

  if(v)
  {
    printf("\nZ_Histo: ");
    for(i=0;i<100;i++)printf("%d ",Z_Histo[i]);
    printf("\n");
  
    printf("MaxHisto:\n");
    for(i=0;i<100;i++)Z_Histo[i]=0;
    for(i=0;i<siglength*5;i++)Z_Histo[(int)MaxCorr[i]]++;
    for(i=0;i<100;i++)printf("%d ",Z_Histo[i]);
    printf("\n");
  }



  //Sorting
  //Reduction to columns
  for(i=0;i<siglength;i++)
  {
    MaxCorr[i]=intmax(MaxCorr[i*5+4],intmax(intmax(MaxCorr[i*5],MaxCorr[i*5+1]),intmax(MaxCorr[i*5+2],MaxCorr[i*5+3])));
  }

  //A histogram of statistical significance
  int *Reihenfolge=Guarded_Malloc(sizeof(int)*siglength*5);
  for(i=0;i<siglength;i++)Reihenfolge[i]=0;
  for(i=0;i<siglength;i++)
  {
    if(MaxCorr[i]<100)Reihenfolge[(int)MaxCorr[i]]++;
    else{Reihenfolge[100]++;}
  }

  if(v)
  {
    printf("The histogram of maximal logprob values for columns:\n");
    fflush(stdout);
    for(i=0;i<101;i++)printf("%d ",Reihenfolge[i]);
    printf("\n");
  }  



  int temporal;
  for(i=0;i<siglength;i++)Reihenfolge[i]=i;

  for(j=0;j<3000;j++)
  {
    for(i=j+1;i<siglength;i++) 
    {
      if(MaxCorr[Reihenfolge[i]]>MaxCorr[Reihenfolge[j]])
      {
        temporal=Reihenfolge[j];
        Reihenfolge[j]=Reihenfolge[i];
        Reihenfolge[i]=temporal;
      }
    }
  }
  if(v)
  {
    printf("These are the strongest Correlations:\n");
    for(j=0;j<1000;j++)printf("%d ",(int)MaxCorr[Reihenfolge[j]]);
    printf("\n");  
  }

  
  //Output of the selection, in the original order, not sorted
  FILE *datei;
  char Chars[6];
  Chars[0]='a';
  Chars[1]='c';
  Chars[2]='g';
  Chars[3]='t';
  Chars[4]='-';
  Chars[5]=' ';
  int *YesNo=Guarded_Malloc(sizeof(int)*siglength);
  int *Choice=Guarded_Malloc(sizeof(int)*1000);


  datei=fopen(outputfile,"w");

  if(NULL == datei){printf("DateiVerbratei!\n");exit(1);}
  //The cho best ones are put into Choice
  for(i=0;i<siglength;i++)YesNo[i]=0;

  while(MaxCorr[Reihenfolge[cho-1]]<c)cho--;
  if(v)printf("chosen: %d\n",cho );

  printf("%d columns were selected as significant.\n",cho );

  for(i=0;i<cho;i++)YesNo[Reihenfolge[i]]=1;
  j=0;  
  for(i=0;i<siglength;i++)
  {
    if(YesNo[i])
    {
      Choice[j]=i;
      j++;
    }
  }  
  //Going through the signatures
  for(j=0;j<signumber;j++)
  {
    for(i=0;i<cho;i++)
    {
      fprintf(datei,"%c", Chars[(int)Signatures[j][Choice[i]]] ); 
    }
    fprintf(datei,"\n");
  }
  fclose(datei);

  return; 
}

void Help()
{
  printf("Usage: ./Correlation MApath\n");
  printf("Flags:\n");
  printf("-b calculates only the correlations of columns with a majority of bases.\n");
  printf("-c <12> changes the correlation cutoff.\n");
  printf("-o file specifies the output file (default SimulatedSignatures).\n");
  printf("-n <300> uses the x best columns above the cutoff to build the signatures.\n");
  printf("-v verbose.\n");
  exit(0);
}




int main(int argc, char *argv[])
{ 

char *MApath_p;
if(argc<2){printf("Usage: ./Correlation MApath\n");exit(0);}
MApath_p=argv[1];

//A few Flags: -b only base columns, -c n with n as corr-cutoff, -o outputfile
int c=12;
int b=0;
int cho=300;
v=0;
char outputfile[]="SimulatedSignatures";
char *output_p=&outputfile[0];
for(int i=1;i<argc;i++)
{
  if(argv[i][0]=='-' && argv[i][1]=='b')
  {
    //printf("%s\n",argv[i]);
    b=1;
  }

  if(argv[i][0]=='-' && argv[i][1]=='c')
  {
    //printf("%s\n",argv[i]);
    if(i+1<argc)c=atoi(argv[i+1]);
    //printf("%d\n",c );
  }

  if(argv[i][0]=='-' && argv[i][1]=='o')
  {
    //printf("%s\n",argv[i]);
    if(i+1<argc)output_p=argv[i+1];
  }

  if(argv[i][0]=='-' && argv[i][1]=='n')
  {
    //printf("%s\n",argv[i]);
    if(i+1<argc)cho=atoi(argv[i+1]);
  }  
  if(argv[i][0]=='-' && argv[i][1]=='v')
  {
    v=1;
  }    
}

if(v)printf("verbose\n");

Einlesen(MApath_p);
if(v)printf("read in\n");

PoissonKorrs(c,b,output_p,cho);

exit(0);

return 0;
}





