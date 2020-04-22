//Secondary function code to run on the GPU with cuda
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>

extern unsigned char *coordinates;

extern unsigned char *elements;

extern size_t nel, nnodes, np, nzmax;

extern char basename;

extern unsigned char *psi;

extern unsigned char *LM;

extern unsigned char *irow;

extern unsigned char *icol;

void num_ElementsNodes(char basename, int myrank)
{
  char fname[30];
  snprintf(fname,100,"%u.vtu",myrank);
  strcat(basename,fname);
  FILE *fp;
  fp = fopen(basename,"r");

  fgets(fp);
  fgets(fp);
  char tline[50];
  fgets(tline,50,fp);
  
  int num = sscanf(tline,"<Piece NumberOfPoints=\"%d\"",int nnodesl);
  int nnodesG = nnodesl;
  
  char str1 = "<Piece NumberOfPoints=\"";
  char str2[30];
  snprintf(str2,100,"%d",nnodesl);
  char str3 = "\"";
  strcat(str1,str2);
  strcat(str1,str3);
  char str4 = " NumberOfCells=\"%d\"";
  strcat(str1,str4);
  int num2 = sscanf(tline,str4,int ncellsl);
  int ncellsG = ncellsl;

  nnodes = nnodesG;
  nel = ncellsG;
}

void read_coordinates(char basename, int myrank, size_t nnodes)
{
  cudaMallocManaged(&coordinates, (nnodes * 3 * sizeof(unsignedchar)));
  
  int count = 0;
  char str = "<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"ascii\">";
  int count = 0;
  printf("got here");
  
  char fname[30];
  snprintf(fname,100,"%u.vtu",myrank);
  strcat(basename,fname);
  FILE *fp;
  fp = fopen(basename,"r");

  char tline[50];
  fgets(tline,50,fp);

  int i;
  int alphabet = 0;
  for (i=0; tline[i]!= '\0'; i++) 
    { 
        // check for alphabets 
        if (isalpha(tline[i]) != 0) 
        {
            alphabet++;
        }
    }
  while(alphabet > 0)
    {
      if(strcmp(tline,str) == 0)
	{
	  break;
	}
      
      fgets(tline,50,fp);

      alphabet = 0;
      for (i=0; tline[i]!= '\0'; i++) 
	{ 
	  // check for alphabets 
	  if (isalpha(tline[i]) != 0) 
	    {
	      alphabet++;
	    }
	}
    }

  for(i=0,i<nnodes,i++)
    {
      fgets(tline,50,fp);
      int num3 = sscanf(tline,"%f %f %f",int x);
      coordinates[count][0] = x[0];
      coordinates[count][1] = x[1];
      coordinates[count][2] = x[2];
      count++;
    }
  fclose(fp);
}

void read_elements(char basename, int myrank, size_t nel)
{
  cudaMallocManaged(&coordinates, (nel * 4 * sizeof(unsignedchar)));
  
  int count = 0;
  char str = "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
  int count = 0;
  printf("got here");
  
  char fname[30];
  snprintf(fname,100,"%u.vtu",myrank);
  strcat(basename,fname);
  FILE *fp;
  fp = fopen(basename,"r");

  char tline[50];
  fgets(tline,50,fp);

  int i;
  int alphabet = 0;
  for (i=0; tline[i]!= '\0'; i++) 
    { 
        // check for alphabets 
        if (isalpha(tline[i]) != 0) 
        {
            alphabet++;
        }
    }
  while(alphabet > 0)
    {
      if(strcmp(tline,str) == 0)
	{
	  break;
	}
      
      fgets(tline,50,fp);

      alphabet = 0;
      for (i=0; tline[i]!= '\0'; i++) 
	{ 
	  // check for alphabets 
	  if (isalpha(tline[i]) != 0) 
	    {
	      alphabet++;
	    }
	}
    }

  for(i=0,i<nnodes,i++)
    {
      fgets(tline,50,fp);
      int num3 = sscanf(tline,"%f %f %f %f",int x);
      coordinates[count][0] = x[0];
      coordinates[count][1] = x[1];
      coordinates[count][2] = x[2];
      coordinates[count][3] = x[3];
      count++;
    }
  fclose(fp);
}

void read_psi(char basename, int myrank, size_t nel)
{
  cudaMallocManaged(&psi, (nel * sizeof(unsignedchar)));
  
  int count = 0;
  char str = "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
  int count = 0;
  printf("got here");
  
  char fname[30];
  snprintf(fname,100,"%u.vtu",myrank);
  strcat(basename,fname);
  FILE *fp;
  fp = fopen(basename,"r");

  char tline[50];
  fgets(tline,50,fp);

  int i;
  int alphabet = 0;
  for (i=0; tline[i]!= '\0'; i++) 
    { 
        // check for alphabets 
        if (isalpha(tline[i]) != 0) 
        {
            alphabet++;
        }
    }
  while(alphabet > 0)
    {
      if(strcmp(tline,str) == 0)
	{
	  break;
	}
      
      fgets(tline,50,fp);

      alphabet = 0;
      for (i=0; tline[i]!= '\0'; i++) 
	{ 
	  // check for alphabets 
	  if (isalpha(tline[i]) != 0) 
	    {
	      alphabet++;
	    }
	}
    }

  for(i=0,i<nnodes,i++)
    {
      fgets(tline,50,fp);
      int x = atoi(tline);
      psi[count] = x;
      count++;
    }
  fclose(fp);
}

struct shapeStruct {
  float N, dN, jac;
};

float struct shapeStruct Struct;

Struct shape(float gp, float xe)
{
  //local coordinate
  float r = gp[0];
  float s = gp[1];
  float t = gp[2];

  //Shape functions
  float N = {r, s, t, 1 - r - s - t};
  float N_r = {1, 0, 0, -1};
  float N_s = {0, 1, 0, -1};
  float N_t = {0 ,0 , 1, -1};

  float x_r = N_r[0]*xe[0,0] + N_r[1]*xe[1,0] + N_r[2]*xe[2,0] + N_r[3]*xe[3,0];
  float x_s = N_s[0]*xe[0,0] + N_s[1]*xe[1,0] + N_s[2]*xe[2,0] + N_s[3]*xe[3,0];
  float x_t = N_t[0]*xe[0,0] + N_t[1]*xe[1,0] + N_t[2]*xe[2,0] + N_t[3]*xe[3,0];

  float y_r = N_r[0]*xe[0,1] + N_r[1]*xe[1,1] + N_r[2]*xe[2,1] + N_r[3]*xe[3,1];
  float y_s = N_s[0]*xe[0,1] + N_s[1]*xe[1,1] + N_s[2]*xe[2,1] + N_s[3]*xe[3,1];
  float y_t = N_t[0]*xe[0,1] + N_t[1]*xe[1,1] + N_t[2]*xe[2,1] + N_t[3]*xe[3,1];
  

__global__ void gen_LMArray(unsigned char elements, size_t nel, int ID)
{
  for(i = blockIdx.x * blockDim.x + threadIdx.x; index < nel; index += (blockDim.x * gridDim.x)))
    {
      for(j=0;j<4;j++)
	{
	  LM[j][i] = ID[elements[i][j]];
	}
    }
}

bool gol_runKernel(unsigned char coordinates, size_t nnodes, float powder_thick, float Tol, unsigned char elements,
		   size_t nel, int ID, ushort threadsCount)
{
  //Get boundary nodes
  int count = 0;
  int i;
  int j;
  int k;
  for(i=0;i<nnodes;i++)
    {
      float z = coordinates[i][2];
      if(fabs(z - powder_thick) < Tol)
	{
	  float fixnodes[count] = i;
	  count++;
	}
    }

  //Assembling ID array
  for(i=0;i<nnodes;i++)
    {
      int ID[i] = 1;
    }
  
  int ndispl = sizeof(fixnodes)/sizeof(fixnodes[0]);
  for(i=0;i<ndispl;i++)
    {
      int nd = fixnodes[i];
      Id[nd] = 0;
    }

  //Fill ID array
  count = 0;
  for(j=0;j<nodes;j++)
    {
      if(ID[j] != 0)
	{
	  count++;
	  ID[j] = count;
	}
    }

  //Generate LM array
  size_t reqBlocksCount = ceil(nel/threadsCount); //number of blocks count for the LM array
  unsigned int blocksCount = (unsigned int)min(65536, (unsigned int)reqBlocksCount); // setting blocks count based on the required blocks.
  gen_LMArray<<<blocksCount, threadsCount>>>(elements, nel, ID);

  int ndof = 0;
  float d;
  for(i = 0;i < nnodes;i++)
    {
      //Displacement Vector
      d[i] = 0;
      if(ID[i] > ndof)
	{
	  ndof = ID[i];
	}
    }
  
  //Compute Sparcity
  nzmax = 0;
  int elem;
  for(elem = 0;elem < nel;elem++)
    {
      for(k = 0;k < 4;k++)
	{
	  int i_index = LM[k][elem];
	  if(i_index > 0)
	    {
	      for(m = 0;m < 4;m++)
		{
		  int j_index = LM[m][elem];
		  if(j_index > 0)
		    {
		      nzmax++;
		    }
		}
	    }
	}
    }

  for(i = 0;i < nzmax;i++)
    {
      irow[i] = 0;
      icol[i] = 0;
    }
  
  count = 0;
  for(elem = 0;elem < nel;elem++)
    {
      for(k = 0;k < 4;k++)
	{
	  i_index = LM[k][elem];
	  if(i_index > 0)
	    {
	      for(m = 0;m < 4;m++)
		{
		  j_index = LM[m][elem];
		  if(j_index > 0)
		    {
		      irow[count] = i_index;
		      icol[count] = j_index;
		      count++;
		    }
		}
	    }
	}
    }

  //Assembling stiffness matrix
  float K;
  float F;
  for(i = 0;i < nzmax;i++)
    {
      K[i] = 0;
    }
  for(i = 0;i < ndof;i++)
    {
      F[i] = 0;
    }
  for(i = 0; i < nel;i++)
    {
      float xe[4][3] = {
			{coordinates[elements[i][0]][0], coordinates[elements[i][0]][1], coordinates[elements[i][0]][2]},
			{coordinates[elements[i][1]][0], coordinates[elements[i][1]][1], coordinates[elements[i][1]][2]},
			{coordinates[elements[i][2]][0], coordinates[elements[i][2]][1], coordinates[elements[i][2]][2]},
			{coordinates[elements[i][3]][0], coordinates[elements[i][3]][1], coordinates[elements[i][3]][2]}
      };
      float Psie = psi[i];
      
  
}
