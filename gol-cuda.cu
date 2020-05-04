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
  extern void num_ElementsNodes(char baseName [80], int myrank, int nnodes);
  
  extern void read_coordinates(char baseName [80], int myrank, int nnodes);

  extern void read_elements(char baseName [80], int myrank, int nel);

  extern void read_psi(char baseName [80], int myrank, int nel);

  extern bool gol_runKernel(float *coordinates, int nnodes, float powderThick,
			    float Tol, float *elements,
			    int nel, int **d_ID, ushort threadsCount,
			    float **d_d, float **d_a_bar);

  extern void gol_freeData();
}

//This is an debugging function that prints out the most recent line the code gets to when running before an error occurs
void printLine(int line)
{
  char fileName[30]; //initializes file buffer
  snprintf(fileName,100,"errorAtLineCuda.txt"); //created file name to be sent to
  FILE *fp;
  fp = fopen(fileName,"w+"); //opens file
  fprintf(fp,"Line is %d.\n",line); //prints the current line number to the file
  fclose(fp); //closes file
}

//This funciton takes in the number of nodes from the user and produces the number of elemements
void num_ElementsNodes(char baseName [80], int myrank, int nnodes)
{
  nel = nnodes * 6; //generates number of elements for specified number of nodes
}

//Generates the dummy coordinate values based off of the specified nodes
void read_coordinates(char baseName [80], int myrank, int nnodes)
{
  int a,b;
  printLine(__LINE__);
  cudaMallocManaged(&coordinates, ((nnodes * 3) * sizeof(float))); //allocates memory
  for (a=0; a<(nnodes); a++)
    {
      for(b=0;b<3;b++)
	{
	  coordinates[a+b*nnodes] = 0; //assigns each coordinate location as zero
	}
    }
}

//Generates dummy element values based off of number of elements
void read_elements(char baseName [80], int myrank, int nel)
{
  cudaMallocManaged(&elements, ((nel * 4) * sizeof(float)));
  int a;
  for (a=0; a<(4*nel); a++)
    {
      elements[a] = 0; //assigns each element value as zero
    }
}

//Generates dummy psi values associated with each element
void read_psi(char baseName [80], int myrank, int nel)
{
  cudaMallocManaged(&psi, (nel * sizeof(float)));
  int a;
  for(a=0;a<nel;a++)
    {
      psi[a] = 0.5; //loads value for psi paramater to each element
    }
}

//Shape function for the FEM approach to deform the mesh
void shape(float gp[3], float xe[12])
{
  int i;
  //local coordinate
  float r = gp[0];
  float s = gp[1];
  float t = gp[2];

  //Shape functions
  N[0] = r;
  N[1] = s;
  N[2] = t;
  N[3] = 1-r-s-t;

  //shape functions at node points
  float N_r[4] = {1, 0, 0, -1};
  float N_s[4] = {0, 1, 0, -1};
  float N_t[4] = {0 ,0 , 1, -1};

  //xyz locations of shape functions at node points
  float x_r = N_r[0]*xe[0+0*4] + N_r[1]*xe[1+0*4] + N_r[2]*xe[2+0*4] + N_r[3]*xe[3+0*4];
  float x_s = N_s[0]*xe[0+0*4] + N_s[1]*xe[1+0*4] + N_s[2]*xe[2+0*4] + N_s[3]*xe[3+0*4];
  float x_t = N_t[0]*xe[0+0*4] + N_t[1]*xe[1+0*4] + N_t[2]*xe[2+0*4] + N_t[3]*xe[3+0*4];

  float y_r = N_r[0]*xe[0+1*4] + N_r[1]*xe[1+1*4] + N_r[2]*xe[2+1*4] + N_r[3]*xe[3+1*4];
  float y_s = N_s[0]*xe[0+1*4] + N_s[1]*xe[1+1*4] + N_s[2]*xe[2+1*4] + N_s[3]*xe[3+1*4];
  float y_t = N_t[0]*xe[0+1*4] + N_t[1]*xe[1+1*4] + N_t[2]*xe[2+1*4] + N_t[3]*xe[3+1*4];

  float z_r = N_r[0]*xe[0+2*4] + N_r[1]*xe[1+2*4] + N_r[2]*xe[2+2*4] + N_r[3]*xe[3+2*4];
  float z_s = N_s[0]*xe[0+2*4] + N_s[1]*xe[1+2*4] + N_s[2]*xe[2+2*4] + N_s[3]*xe[3+2*4];
  float z_t = N_t[0]*xe[0+2*4] + N_t[1]*xe[1+2*4] + N_t[2]*xe[2+2*4] + N_t[3]*xe[3+2*4];

  //determinant of the jacobian
  float jacDet = x_r*(y_s*z_t - y_t*z_s) - x_s*(y_r*z_t - y_t*z_r) + x_t*(y_r*z_s - y_s*z_r);

  //jacobian
  jac = abs(jacDet);

  //inverse of the jacobian
  float inv_jac[3][3] = {(y_s*z_t - y_t*z_s)/jacDet, (x_t*z_s - x_s*z_t)/jacDet, (x_s*y_t - x_t*y_s)/jacDet,
                         (y_t*z_r - y_r*z_t)/jacDet, (x_r*z_t - x_t*z_r)/jacDet, (x_t*y_r - x_r*y_t)/jacDet,
                         (y_r*z_s - y_s*z_r)/jacDet, (x_s*z_r - x_r*z_s)/jacDet, (x_r*y_s - x_s*y_r)/jacDet};

  //derivative of the shape function
  for(i=0;i<4;i++)
    {
      dN[i+0*4] = N_r[i]*inv_jac[0][0] + N_s[i]*inv_jac[1][0] + N_t[i]*inv_jac[2][0];
      dN[i+1*4] = N_r[i]*inv_jac[0][1] + N_s[i]*inv_jac[1][1] + N_t[i]*inv_jac[2][1];
      dN[i+2*4] = N_r[i]*inv_jac[0][2] + N_s[i]*inv_jac[1][2] + N_t[i]*inv_jac[2][2];
    }
}

void weakform(float xe[12], float Psie, float porosity)
{
  int i;
  int j;
  int k;
  int l;

  // 1 point formula - degree of precision 1
  float gp[3] = {0.25, 0.25, 0.25};

  //weights
  int w = 1;

  int ngp = 1;

  for(i=0;i<16;i++)
    {
      ke[i] = 0; //initialize ke values
    }
  
  for(i=0;i<4;i++)
    {
      fe[i] = 0; //initialize fe values
    }

  //stress strain displacement matrix
  float B[4] = {0,0,0,0};
  
  //loop over gauss points
  shape(gp,xe); //use shape function outlined above

  for(i=0;i<ngp;i++)
    {
      float por = porosity;

      //current location in z of the point
      float z = N[0]*xe[0+2*4] + N[1]*xe[1+2*4] + N[2]*xe[2+2*4] + N[3]*xe[3+2*4];
      if( z < powderThick)
        {
          por = 0.0; //porosit of the substrate is zero since it's a solid
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
}

//function that conversts the solution of the system of equations to the displacements of the powder layer node points.
__global__ void gen_corrector(int nnodes, int *d_ID, float *d_a_bar, float powderThick, float *d_d)
{
  int i, index;

  //Use CUDA to run the displacment converison calculation
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
        }
    }
}

//runs the majority of the displacment calculation and the CUDA kernel outlined above
bool gol_runKernel(float *coordinates, int nnodes, float powderThick,float Tol, float *elements, int nel, int **d_ID, ushort threadsCount, float **d_d, float **d_a_bar)
{
  printLine(__LINE__);
  //Get boundary nodes
  int count = 0;
  int i;
  int j;
  int *fixnodes;
  float z;

  //set fixnodes for boundary nodes to close to the boundary of the powder to the subsrate
  for(i=0;i<nnodes;i++)
    {
      z = coordinates[i+2*nnodes];
      if(fabs(z - powderThick) < Tol)
        {
          fixnodes[count] = i;
          count++;
        }
    }
  printLine(__LINE__);
  cudaMallocManaged(&ID, (nnodes * sizeof(int))); //allocate space
  cudaMallocManaged(&d, (nnodes * sizeof(float))); //allocate space
  printLine(__LINE__);
  
  //Assembling ID array
  float ID[nnodes];
  for(i=0;i<nnodes;i++)
    {
      ID[i] = 1;
    }
  printLine(__LINE__);
  int ndispl = sizeof(fixnodes)/sizeof(fixnodes[0]);

  int nd;
  int g;
  printf("ndispl is %d\n",ndispl);
  printLine(__LINE__);

  //Adjust ID array based on the nd values
  for(g=0; g<ndispl; g++)
    {
      nd = g;
      ID[nd] = 0;
    }
  printLine(__LINE__);
  
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
  printLine(__LINE__);
  int ndof = 1;

  //initialize values for the displacement vector
  for(i = 0;i < nnodes;i++)
    {
      //Displacement Vector
      d[i] = 0;
      if(ID[i] > ndof)
        {
          ndof = ID[i];
        }
    }
  printf("ndof is %d\n",ndof);
  printLine(__LINE__);

  //For future work implement solver for solving system of equation M * a_bar = F
  cudaMallocManaged(&a_bar, (ndof * sizeof(float))); //allocate space
  printLine(__LINE__);
  
  float F[ndof];
  for(i = 0;i < ndof;i++)
    {
      F[i] = 0; //b value to system of equations
    }
  printLine(__LINE__);
  float a_bar[ndof];

  for(i=0; i<ndof;i++)
  {
    a_bar[i] = F[i]; //set a_bar value
  }
  printLine(__LINE__);
  
  //Change small values to zero
  for(i=0;i<ndof;i++)
    {
      if(a_bar[i] < 0.0000001)
        {
          a_bar[i] = 0;
        }
    }
  printLine(__LINE__);
  
  //Corrector phase to be done in cuda 
  size_t reqBlocksCount2 = ceil(nnodes/threadsCount); //number of blocks count for the LM array
  unsigned int blocksCount2 = (unsigned int)min(65536, (unsigned int)reqBlocksCount2); //Blocks count based on required blocks
  gen_corrector<<<blocksCount2, threadsCount>>>(nnodes, *d_ID, *d_a_bar, powderThick, *d_d);
  cudaDeviceSynchronize();
  printLine(__LINE__);
  for(i=0;i<nnodes;i++)
    {
      coordinates[i+2*nnodes] = coordinates[i+2*nnodes] + d[i]; //deform coordinates based on d
    }
  printLine(__LINE__);
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
