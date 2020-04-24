//Secondary function code to run on the GPU with cuda
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include<unistd.h>
#include<stdbool.h>
#include<cuda.h>
#include<cuda_runtime.h>

extern unsigned char *coordinates;

extern unsigned char *elements;

extern size_t nel, nnodes, np, nzmax;

extern char basename;

extern unsigned char *psi;

extern unsigned char *LM;

extern unsigned char *irow;

extern unsigned char *icol;

extern "C"
{
  extern void num_ElementsNodes(char basename, int myrank);

  extern void read_coordinates(char basename, int myrank, size_t nnodes);

  extern void read_elements(char basename, int myrank, size_t nel);

  extern void read_psi(char basename, int myrank, size_t nel);

  extern bool gol_runKernel(unsigned char coordinates, size_t nnodes, float powder_thick, float Tol, unsigned char elements,
			  size_t nel, ushort threadsCount);
}

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
  //cudaMallocManaged(&coordinates, (nnodes * 3 * sizeof(unsignedchar)));
  
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
  //cudaMallocManaged(&coordinates, (nel * 4 * sizeof(unsignedchar)));
  
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
  //cudaMallocManaged(&psi, (nel * sizeof(unsignedchar)));
  
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
  Struct s;
  int i;
  //local coordinate
  float r = gp[0];
  float s = gp[1];
  float t = gp[2];

  //Shape functions
  float N = {r, s, t, 1 - r - s - t};
  float N_r = {1, 0, 0, -1};
  float N_s = {0, 1, 0, -1};
  float N_t = {0 ,0 , 1, -1};

  float x_r = N_r[0]*xe[0][0] + N_r[1]*xe[1][0] + N_r[2]*xe[2][0] + N_r[3]*xe[3][0];
  float x_s = N_s[0]*xe[0][0] + N_s[1]*xe[1][0] + N_s[2]*xe[2][0] + N_s[3]*xe[3][0];
  float x_t = N_t[0]*xe[0][0] + N_t[1]*xe[1][0] + N_t[2]*xe[2][0] + N_t[3]*xe[3][0];

  float y_r = N_r[0]*xe[0][1] + N_r[1]*xe[1][1] + N_r[2]*xe[2][1] + N_r[3]*xe[3][1];
  float y_s = N_s[0]*xe[0][1] + N_s[1]*xe[1][1] + N_s[2]*xe[2][1] + N_s[3]*xe[3][1];
  float y_t = N_t[0]*xe[0][1] + N_t[1]*xe[1][1] + N_t[2]*xe[2][1] + N_t[3]*xe[3][1];

  float z_r = N_r[0]*xe[0][2] + N_r[1]*xe[1][2] + N_r[2]*xe[2][2] + N_r[3]*xe[3][2];
  float z_s = N_s[0]*xe[0][2] + N_s[1]*xe[1][2] + N_s[2]*xe[2][2] + N_s[3]*xe[3][2];
  float z_t = N_t[0]*xe[0][2] + N_t[1]*xe[1][2] + N_t[2]*xe[2][2] + N_t[3]*xe[3][2];

  float jacobian[3][3] = {x_r,x_s,x_t,y_r,y_s,y_t,z_r,z_s,z_t};
  float jacDet = x_r*(y_s*z_t - y_t*z_s) - x_s*(y_r*z_t - t_t*z_r) + x_t*(y_r*z_s - y_s*z_r);
  float jac = abs(jacDet);

  //Check Jacobian
  if(jac <= 0.0)
    {
      fprintf(strerr,"Negative jacobian, element too distorted!\n");
    }
  //Take the inverse of the Jacobian
  float inv_jac[3][3] = {(y_s*z_t - y_t*z_s)/jacDet, (x_t*z_s - x_s*z_t)/jacDet, (x_s*y_t - x_t*y_s)/jacDet,
			 (y_t*z_r - y_r*z_t)/jacDet, (x_r*z_t - x_t*z_r)/jacDet, (x_t*y_r - x_r*y_t)/jacDet,
			 (y_r*z_s - y_s*z_r)/jacDet, (x_s*z_r - x_r*z_s)/jacDet, (x_r*y_s - x_s*y_r)/jacDet};
  float dN;
  for(i=0;i<4,i++)
    {
      dN[i][0] = N_r[i]*inv_jac[0][0] + N_s[i]*inv_jac[1][0] + N_t[i]*inv_jac[2][0];
      dN[i][1] = N_r[i]*inv_jac[0][1] + N_s[i]*inv_jac[1][1] + N_t[i]*inv_jac[2][1];
      dN[i][2] = N_r[i]*inv_jac[0][2] + N_s[i]*inv_jac[1][2] + N_t[i]*inv_jac[2][2];
    }
  s.N = N;
  s.dN = dN;
  s.jac = jac;

  return s;
}

struct weakformStruct {
  float ke, fe;
};

float weakformStruct Struct;

Struct weakform(float xe, float Psie, float porosity)
{
  int i, j, k;
  // 1 point formula - degree of precision 1
  float gp = {0.25, 0.25, 0.25};
  int w = 1;

  int ngp = 1;

  //initialize stiffness matrix
  float ke[4][4] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  //right hand size
  float fe[4][1] = {0,0,0,0};

  //stress strain displacement matrix
  float B[1][4] = {0,0,0,0};
  //loop over gauss points
  Struct result = shape(gp,xe);
  float N = result.N;
  float dN = result.dN;
  float jac = result.jac;
  for(i=0;i<ngp;i++)
    {
      float por = porosity;
      float z = N * {xe[0][2], xe[1][2], xe[2][2], xe[3][2]};
      if( z < powderThick)
	{
	  por = 0.0;
	}
      for(j=0;j<4;j++)
	{
	  B[j] = dN[j][3];
	}
      //Transpose of N
      float Ntr[4][1] = {N[0],N[1],N[2],N[3]};
      //fill k
      ke = ke + Ntr * B * w[i] * jac;
      //fill fe
      fe = fe - Ntr * ((por * Psie)/(1 - por * (1 - Psie)))*w[i]*jac;
    }
  s.ke = ke;
  s.fe = fe;

  return s;
}


__global__ void gen_LMArray(unsigned char elements, size_t nel, int ID)
{
  for(i = blockIdx.x * blockDim.x + threadIdx.x; i < nel; i += (blockDim.x * gridDim.x))
    {
      for(j=0;j<4;j++)
	{
	  LM[j][i] = ID[elements[i][j]];
	}
    }
}

__global__ float gen_corrector(size_t nnodes, int ID, float a_bar, float powderThick)
{
  int i, index;
  float d;
  for(i = blockIdx.x * blockDim.x + threadIdx.x; i < nnodes; i += (blockDim.x * gridDim.x))
    {
      index = ID[i];
      if(index != 0)
	{
	  d[i] = a_bar[index];
	  if(d[i] > powderThick)
	    {
	      d[i] = 0;
	    }
	  if(isdigit(round(d[i])) == 0)
	    {
	      d[i] = 0;
	    }
	}
    }

  return d
}

bool gol_runKernel(unsigned char coordinates, size_t nnodes, float powder_thick, float Tol, unsigned char elements,
		   size_t nel, ushort threadsCount)
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

  cudaMallocManaged(&ID, (nnodes * sizeof(int)));
  
  //Assembling ID array
  int ID;
  for(i=0;i<nnodes;i++)
    {
      ID[i] = 1;
    }
  
  int ndispl = sizeof(fixnodes)/sizeof(fixnodes[0]);
  int nd;
  for(i=0;i<ndispl;i++)
    {
      nd = fixnodes[i];
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
  
  cudaMallocManaged(&LM, (nel * sizeof(unsigned char)));
  
  //Generate LM array
  size_t reqBlocksCount = ceil(nel/threadsCount); //number of blocks count for the LM array
  unsigned int blocksCount = (unsigned int)min(65536, (unsigned int)reqBlocksCount); // setting blocks count based on the required blocks.
  gen_LMArray<<<blocksCount, threadsCount>>>(elements, nel, ID);
  cudaDeviceSynchronize();
  
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
  int elem, i_index, j_index;
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
  float xe;
  float Psie;
  Struct result;
  float ke;
  float fe;
  count = 0;
  for(i = 0;i < nzmax;i++)
    {
      K[i] = 0;
    }
  for(i = 0;i < ndof;i++)
    {
      F[i][0] = 0;
    }
  for(i = 0; i < nel;i++)
    {
      xe[4][3] = {
		  {coordinates[elements[i][0]][0], coordinates[elements[i][0]][1], coordinates[elements[i][0]][2]},
		  {coordinates[elements[i][1]][0], coordinates[elements[i][1]][1], coordinates[elements[i][1]][2]},
		  {coordinates[elements[i][2]][0], coordinates[elements[i][2]][1], coordinates[elements[i][2]][2]},
		  {coordinates[elements[i][3]][0], coordinates[elements[i][3]][1], coordinates[elements[i][3]][2]}
      };
      Psie = psi[i];
      result = weakform(xe,Psie,porosity);
      ke = result.ke;
      fe = result.fe;
      for(j=0;j<4;j++)
	{
	  i_index = LM[j][i];
	  if(i_index > 0)
	    {
	      F[i_index] = F[i_index] + fe[j];
	      for(k=0;k<4;k++)
		{
		  j_index = LM[k][i];
		  if(j_index > 0)
		    {
		      K[count] = K[count] + ke[j][k];
		    }
		}
	    }
	}
    }

  float M;
  //Ensamble sparse matrix
  for(i=0;i<ndof;i++)
    {
      for(j=0;j<ndof;j++)
	{
	  M[i][j] = 0;
	}
    }
  for(i=0;i<nzmax;i++)
    {
      M[irow[i]][icol[i]] = M[irow[i]][icol[i]] + K[i];
    }
  //For future work implement solver for solving system of equation M * a_bar = F
  cudaMallocManaged(&a_bar, (ndof * sizeof(float)));
  
  float a_bar;
  a_bar = F;

  //Change small values to zero
  for(i=0;i<ndof;i++)
    {
      if(a_bar[i] < 0.0000001)
	{
	  a_bar[i] = 0;
	}
    }

  //Corrector phase to be done in cuda
  reqBlocksCount = ceil(nnodes/threadsCount); //number of blocks count for the LM array
  blocksCount = (unsigned int)min(65536, (unsigned int)reqBlocksCount); // setting blocks count based on the required blocks.
  d = gen_corrector<<<blocksCount, threadsCount>>>(nnodes, ID, a_bar, powderThick);
  cudaDeviceSynchronize();
  
  for(i=0;i<nnodes;i++)
    {
      coordinates[i][2] = coordinates[i][2] + d[i];
    }
}
