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

float nnodesl;

extern int *ID;

extern float *coordinates;

extern float *elements;

extern int nel, nnodes, np, nzmax;

extern char baseName[80];

extern float *psi;

extern float *LM;

extern int *irow;

extern int *icol;

extern float powderThick;

float *N;

float *dN;

float jac;

float *ke;

float *fe;

extern float *d;

extern float *a_bar;

extern "C"
{
  extern void num_ElementsNodes(char baseName [80], int myrank);

  extern int offsetCalc(char baseName [80], int numranks, int myrank);

  extern void read_coordinates(char baseName [80], int myrank, int nnodes);

  extern void read_elements(char baseName [80], int myrank, int nel);

  extern void read_psi(char baseName [80], int myrank, int nel);

  extern bool gol_runKernel(float *coordinates, int nnodes, float powderThick,
			    float Tol, float *elements,
			    int nel, int **d_ID, ushort threadsCount,
			    float **d_d, float **d_a_bar);

  extern void gol_freeData();
}

void printLine(int line)
{
  char fileName[30];
  snprintf(fileName,100,"errorAtLineCuda.txt");
  FILE *fp;
  fp = fopen(fileName,"w+");
  fprintf(fp,"Line is %d.\n",line);
  fclose(fp);
}

void num_ElementsNodes(char baseName [80], int myrank)
{
  printLine(__LINE__);
  char fname [100];
  snprintf(fname,100,"%u.vtu",myrank);
  strcat(baseName, fname);
  FILE *fp;
  fp = fopen(baseName,"r");
  printLine(__LINE__);
  
  fgets(baseName, 47, fp);
  char tline[50];
  fgets(tline, 50, fp);
  printLine(__LINE__);
  
  int num = sscanf(tline,"<Piece NumberOfPoints=\"%d\"",&nnodesl);
  int nnodesG = nnodesl;
  printLine(__LINE__);
  
  char str1 [24] = "<Piece NumberOfPoints=\"";
  char str2 [30];
  snprintf(str2,100,"%d",nnodesl);
  char str3 [3] = "\"";
  strcat(str1,str2);
  strcat(str1,str3);
  char str4 [100] = "NumberOfCells=\"%d\"";
  strcat(str1,str4);
  int ncellsl;
  int num2 = sscanf(tline,str4, &ncellsl);
  int ncellsG = ncellsl;
  printLine(__LINE__);
  
  nnodes = nnodesG;
  nel = ncellsG;
}

int offsetCalc(char baseName [80], int numranks, int myrank)
{
  int i;
  char fname[30];
  FILE *fp;
  char tline[50];
  float nnodesG[numranks];
  float offset[numranks];
  for(i=0;i<numranks;i++)
    {
      snprintf(fname,100,"%u.vtu",myrank);
      strcat(baseName,fname);
      fp = fopen(baseName,"r");

      fgets(tline,50,fp);

      float num = sscanf(tline,"<Piece NumberOfPoints=\"%d\"", nnodesl);
      nnodesG[i] = nnodesl;
    }

  for(i=0;i<numranks;i++)
    {
      offset[i] = 0;
    }

  for(i=1;i<numranks;i++)
    {
      offset[i] = offset[i-1] + nnodesG[i-1];
    }
    return offset[numranks];
}

void read_coordinates(char baseName [80], int myrank, int nnodes)
{
  int a,b;
  //float coordinates[nnodes][3];
  cudaMallocManaged(&coordinates, ((nnodes * 3) * sizeof(float)));
  //size_t pitch;
  //cudaMallocPitch(&coordinates, &pitch, sizeof(float)*3, nnodes);
  for (a=0; a<(nnodes); a++)
    {
      for(b=0;b<3;b++)
	{
	  coordinates[a+b*nnodes] = 0;
	}
    }


  int count = 0;
  char str[100] = "<DataArray type=\"Float64\" Name=\"coordinates\" NumberOfComponents=\"3\" format=\"ascii\">";
  printf("got here");

  char fname[30];
  snprintf(fname,100,"%u.vtu",myrank);
  strcat(baseName,fname);
  FILE *fp;
  fp = fopen(baseName,"r");

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
  //coordinates[nnodes][3];
  for(i=0;i<nnodes;i++)
    {
      float x,y,z;
      //float coordinates[nnodes][3];
      fgets(tline,50,fp);
      int num3 = sscanf(tline,"%f %f %f", &x,&y,&z);
      coordinates[count] = x;
      coordinates[count+nnodes] = y;
      coordinates[count+2*nnodes] = z;
      count++;
    }
  fclose(fp);
}

void read_elements(char baseName [80], int myrank, int nel)
{
  cudaMallocManaged(&elements, ((nel * 4) * sizeof(float)));
  int a;
  for (a=0; a<(4*nel); a++)
    {
      elements[a] = 0;
    }

  int count = 0;
  char str[80] = "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
  printf("got here");

  char fname[30];
  snprintf(fname,100,"%u.vtu",myrank);
  strcat(baseName,fname);
  FILE *fp;
  fp = fopen(baseName,"r");

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

  for(i=0;i<nel;i++)
    {
      fgets(tline,50,fp);
      float x,y,z,w;
      //float elements[nel][4];
      int num3 = sscanf(tline,"%f %f %f %f", &x,&y,&z,&w);
      elements[count] = x;
      elements[count+nel] = y;
      elements[count+2*nel] = z;
      elements[count+3*nel] = w;
      count++;
    }
  fclose(fp);
}

void read_psi(char baseName [80], int myrank, int nel)
{
  cudaMallocManaged(&psi, (nel * sizeof(float)));

  int count = 0;
  char str[100] = "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">";
  printf("got here");

  char fname[30];
  snprintf(fname,100,"%u.vtu",myrank);
  strcat(baseName,fname);
  FILE *fp;
  fp = fopen(baseName,"r");

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

  for(i=0;i<nel;i++)
    {
      fgets(tline,50,fp);
      float x = atoi(tline);
      psi[count] = x;
      count++;
    }
  fclose(fp);
}

//struct shapeStruct {
//  float NS, dNS, jacS;
//};

//typedef struct shapeStruct Struct;

//Struct shape(float gp[3], float xe[4][3])
//void shape(float gp[3], float xe[4][3], float *add_N[], float *add_dN[][], float *add_jac)
void shape(float gp[3], float xe[12])
{
  //struct s;
  int i;
  //float xe[4][3];
  //float gp[3];
  //local coordinate
  float r = gp[0];
  float s = gp[1];
  float t = gp[2];

  //Shape functions
  //N[4] = {r, s, t, 1 - r - s - t};
  N[0] = r;
  N[1] = s;
  N[2] = t;
  N[3] = 1-r-s-t;
  float N_r[4] = {1, 0, 0, -1};
  float N_s[4] = {0, 1, 0, -1};
  float N_t[4] = {0 ,0 , 1, -1};

  float x_r = N_r[0]*xe[0+0*4] + N_r[1]*xe[1+0*4] + N_r[2]*xe[2+0*4] + N_r[3]*xe[3+0*4];
  float x_s = N_s[0]*xe[0+0*4] + N_s[1]*xe[1+0*4] + N_s[2]*xe[2+0*4] + N_s[3]*xe[3+0*4];
  float x_t = N_t[0]*xe[0+0*4] + N_t[1]*xe[1+0*4] + N_t[2]*xe[2+0*4] + N_t[3]*xe[3+0*4];

  float y_r = N_r[0]*xe[0+1*4] + N_r[1]*xe[1+1*4] + N_r[2]*xe[2+1*4] + N_r[3]*xe[3+1*4];
  float y_s = N_s[0]*xe[0+1*4] + N_s[1]*xe[1+1*4] + N_s[2]*xe[2+1*4] + N_s[3]*xe[3+1*4];
  float y_t = N_t[0]*xe[0+1*4] + N_t[1]*xe[1+1*4] + N_t[2]*xe[2+1*4] + N_t[3]*xe[3+1*4];

  float z_r = N_r[0]*xe[0+2*4] + N_r[1]*xe[1+2*4] + N_r[2]*xe[2+2*4] + N_r[3]*xe[3+2*4];
  float z_s = N_s[0]*xe[0+2*4] + N_s[1]*xe[1+2*4] + N_s[2]*xe[2+2*4] + N_s[3]*xe[3+2*4];
  float z_t = N_t[0]*xe[0+2*4] + N_t[1]*xe[1+2*4] + N_t[2]*xe[2+2*4] + N_t[3]*xe[3+2*4];

  float jacDet = x_r*(y_s*z_t - y_t*z_s) - x_s*(y_r*z_t - y_t*z_r) + x_t*(y_r*z_s - y_s*z_r);
  jac = abs(jacDet);

  //Check Jacobian
  //if(jac =< 0.0)
  //  {
  //    fprintf(stderr, "Negative jacobian, element too distorted!\n");
  //  }
  //Take the inverse of the Jacobian
  float inv_jac[3][3] = {(y_s*z_t - y_t*z_s)/jacDet, (x_t*z_s - x_s*z_t)/jacDet, (x_s*y_t - x_t*y_s)/jacDet,
                         (y_t*z_r - y_r*z_t)/jacDet, (x_r*z_t - x_t*z_r)/jacDet, (x_t*y_r - x_r*y_t)/jacDet,
                         (y_r*z_s - y_s*z_r)/jacDet, (x_s*z_r - x_r*z_s)/jacDet, (x_r*y_s - x_s*y_r)/jacDet};
  //dN[4][3];
  for(i=0;i<4;i++)
    {
      dN[i+0*4] = N_r[i]*inv_jac[0][0] + N_s[i]*inv_jac[1][0] + N_t[i]*inv_jac[2][0];
      dN[i+1*4] = N_r[i]*inv_jac[0][1] + N_s[i]*inv_jac[1][1] + N_t[i]*inv_jac[2][1];
      dN[i+2*4] = N_r[i]*inv_jac[0][2] + N_s[i]*inv_jac[1][2] + N_t[i]*inv_jac[2][2];
    }
  //s.NS = N;
  //s.dNS = dN;
  //s.jacS = jac;

  //*add_N = N[];
  //*add_dN = dN[][];
  //*add_jac = jac;

  //return s;
}

//struct weakformStruct {
//  float ke, fe;
//};

//typedef struct weakformStruct Struct1;

//Struct1 weakform(float xe[4][3], float Psie, float porosity);
void weakform(float xe[12], float Psie, float porosity)
{
  int i;
  int j;
  int k;
  int l;
  //float xe[4][3];
  //float N, dN, jac;
  // 1 point formula - degree of precision 1
  float gp[3] = {0.25, 0.25, 0.25};
  int w = 1;

  int ngp = 1;

  //initialize stiffness matrix
  //float ke[4][4] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  for(i=0;i<16;i++)
    {
      ke[i] = 0;
    }
  
  //right hand size
  //float fe[4][1] = {0,0,0,0};
  for(i=0;i<4;i++)
    {
      fe[i] = 0;
    }

  //stress strain displacement matrix
  float B[4] = {0,0,0,0};
  //loop over gauss points
  //struct result = shape(gp, xe[4][3]);
  //shape(gp,xe,&N,&dN,&jac);
  shape(gp,xe);
  //float N = result.N;
  //float dN = result.dN;
  //float jac = result.jac;
  for(i=0;i<ngp;i++)
    {
      float por = porosity;
      //float z = N * {xe[0][2], xe[1][2], xe[2][2], xe[3][2]};
      float z = N[0]*xe[0+2*4] + N[1]*xe[1+2*4] + N[2]*xe[2+2*4] + N[3]*xe[3+2*4];
      if( z < powderThick)
        {
          por = 0.0;
        }
      for(j=0;j<4;j++)
        {
          B[j] = dN[j*2*4];
        }
      //Transpose of N
      float Ntr[4] = {N[0],N[1],N[2],N[3]};
      //fill k
      for(k=0;k<4;k++)
	{
	  for(l=0;l<4;l++)
	    {
	      ke[k+1] = ke[k+l] + Ntr[l] * B[k] * w * jac;
	    }
	}
      //fill fe
      for(k=0;k<4;k++)
	{
	  fe[k] = fe[k] - Ntr[k] * ((por * Psie)/(1 - por * (1 - Psie)))*w*jac;
	}
    }
  //float s.ke = ke;
  //float s.fe = fe;

  //return s;
}

/*
__global__ void gen_LMArray(float elements, int nel)
{
  int i,j;
  for(i = blockIdx.x * blockDim.x + threadIdx.x; i < nel; i += (blockDim.x * gridDim.x))
    {
      for(j=0;j<4;j++)
        {
          LM[j+i*nel] = ID[elements[i+j*nel]];
        }
    }
}
*/

__global__ void gen_corrector(int nnodes, int *d_ID, float *d_a_bar, float powderThick, float *d_d)
{
  int i, index;
  //float d;
  for(i = blockIdx.x * blockDim.x + threadIdx.x; i < nnodes; i += (blockDim.x * gridDim.x))
    {
      index = d_ID[i];
      if(index != 0)
        {
          d_d[i] = d_a_bar[index];
          if(d_d[i] > powderThick)
            {
              d_d[i] = 0;
            }
          //if(isdigit(round(d_d[i])) == 0)
          //  {
          //    d_d[i] = 0;
          //  }
        }
    }

  //return d
}

bool gol_runKernel(float *coordinates, int nnodes, float powderThick,float Tol, float *elements, int nel, int **d_ID, ushort threadsCount, float **d_d, float **d_a_bar)
{
  //Get boundary nodes
  int count = 0;
  int i;
  int j;
  //int k;
  //int m;
  int *fixnodes;
  //float *xe;
  for(i=0;i<nnodes;i++)
    {
      float z = coordinates[i+2*nnodes];
      if(fabs(z - powderThick) < Tol)
        {
          fixnodes[count] = i;
          count++;
        }
    }

  cudaMallocManaged(&ID, (nnodes * sizeof(int)));
  cudaMallocManaged(&d, (nnodes * sizeof(float)));

  //Assembling ID array
  float ID[nnodes];
  for(i=0;i<nnodes;i++)
    {
      ID[i] = 1;
    }

  int ndispl = sizeof(fixnodes)/sizeof(fixnodes[0]);
  int nd;
  int g;
  for(g=0; g<ndispl; g++)
    {
      nd = fixnodes[g];
      ID[nd] = 0;
    }

  //Fill ID array
  count = 0;
  for(j=0;j<nnodes;j++)
    {
      if(ID[j] != 0)
        {
          count++;
          ID[j] = count;
        }
    }

  int ndof = 0;
  //float d;
  for(i = 0;i < nnodes;i++)
    {
      //Displacement Vector
      d[i] = 0;
      if(ID[i] > ndof)
        {
          ndof = ID[i];
        }
    }
  
  /*
  cudaMallocManaged(&LM, (nel * sizeof(unsigned char)));

  //Generate LM array
  size_t reqBlocksCount = ceil(nel/threadsCount); //number of blocks count for the LM array
  unsigned int blocksCount = (unsigned int)min(65536, (unsigned int)reqBlocksCount); // setting blocks count based on the required blocks.
  int gen_LMArray<<<blocksCount, threadsCount>>>(elements, nel, ID);
  cudaDeviceSynchronize();

  //Compute Sparcity
  nzmax = 0;
  int elem, i_index, j_index;
  for(elem = 0;elem < nel;elem++)
    {
      for(k = 0;k < 4;k++)
        {
          i_index = LM[k+elem*nel];
          if(i_index > 0)
            {
              for(m = 0;m < 4;m++)
                {
                  j_index = LM[m+elem*nel];
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
  for(elem = 0; elem < nel; elem++)
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
  float F[ndof];
  float Psie;
  struct result;
  count = 0;
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
      for(j=0;j<3;j++)
	{
	  //xe[4][3] = {
	  //	      {coordinates[elements[i][0]][0], coordinates[elements[i][0]][1], coordinates[elements[i][0]][2]},
	  //	      {coordinates[elements[i][1]][0], coordinates[elements[i][1]][1], coordinates[elements[i][1]][2]},
	  //	      {coordinates[elements[i][2]][0], coordinates[elements[i][2]][1], coordinates[elements[i][2]][2]},
	  //	      {coordinates[elements[i][3]][0], coordinates[elements[i][3]][1], coordinates[elements[i][3]][2]}
	  //};
	  xe[0+j*4] = coordinates[elements[i+0*nel]+j*nnodes];
	  xe[1+j*4] = coordinates[elements[i+1*nel]+j*nnodes];
	  xe[2+j*4] = coordinates[elements[i+2*nel]+j*nnodes];
	  xe[3+j*4] = coordinates[elements[i+3*nel]+j*nnodes];
	}
      Psie = psi[i];
      result = weakform(xe,Psie,porosity);
      //ke = result.ke;
      //fe = result.fe;
      for(j=0;j<4;j++)
        {
          i_index = LM[j+i*nel];
          if(i_index > 0)
            {
              F[i_index] = F[i_index] + fe[j];
              for(k=0;k<4;k++)
                {
                  j_index = LM[k][i];
                  if(j_index > 0)
                    {
                      K[count] = K[count] + ke[j+k*4];
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
          M[i+j*ndof] = 0;
        }
    }
  for(i=0;i<nzmax;i++)
    {
      M[irow[i]+icol[i]*nzmax] = M[irow[i]+icol[i]*nzmax] + K[i];
    }

  */
  //For future work implement solver for solving system of equation M * a_bar = F
  cudaMallocManaged(&a_bar, (ndof * sizeof(float)));

  float *F;
  for(i = 0;i < ndof;i++)
    {
      F[i] = 0;
    }
  
  float *a_bar;

  for(i=0; i<ndof;i++)
  {
      a_bar[i] = F[i];
  }

  //Change small values to zero
  for(i=0;i<ndof;i++)
    {
      if(a_bar[i] < 0.0000001)
        {
          a_bar[i] = 0;
        }
    }

  //Corrector phase to be done in cuda 
  size_t reqBlocksCount2 = ceil(nnodes/threadsCount); //number of blocks count for the LM array
  unsigned int blocksCount2 = (unsigned int)min(65536, (unsigned int)reqBlocksCount2);
  gen_corrector<<<blocksCount2, threadsCount>>>(nnodes, *d_ID, *d_a_bar, powderThick, *d_d);
  cudaDeviceSynchronize();

  for(i=0;i<nnodes;i++)
    {
      coordinates[i+2*nnodes] = coordinates[i+2*nnodes] + d[i];
    }

  return 0;
}

//free the cuda data that is allocated.
void gol_freeData()
{
 cudaFree(coordinates);
 cudaFree(elements);
 cudaFree(psi);
 cudaFree(d);
 cudaFree(ID);
 cudaFree(a_bar);
}
