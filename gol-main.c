//Main function code to run on the CPU
//Jacob Pessin & Alex Vest, April 2020 Parallel Programming final project
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include<unistd.h>
#include<stdbool.h>
#include<mpi.h>

int *ID;

float *d;

float *a_bar;

float *coordinates;

float *elements;

int nel, nnodes, nzmax;

char baseName[80]  = "/Powder-Layer-Consolidation/Layer_270_09_10/0/";

float *psi;

float *LM;

int *irow;

int *icol;

float porosity = 0.8;

float powderThick = 28.5;

float Tol = 2.5;

size_t np = 4;

//extern functions from the gol-cuda file
extern void num_ElementsNodes(char baseName[80], int myrank);

extern int offsetCalc(char baseName[80], int numranks, int myrank);

extern void read_coordinates(char baseName[80], int myrank, int nnodes);

extern void read_elements(char baseName[80], int myrank, int nel);

extern void read_psi(char baseName[80], int myrank, int nel);

extern bool gol_runKernel(float *coordinates, int nnodes, float powderThick,
			  float Tol, float *elements,
                          int nel, int **d_ID, ushort threadsCount,
			  float **d_d, float **d_a_bar);

extern void gol_freeData();

void printLine(int line)
{
  char fileName[30];
  snprintf(fileName,100,"errorAtLine.txt");
  FILE *fp;
  fp = fopen(fileName,"w+");
  fprintf(fp,"Line is %d.\n",line);
  fclose(fp);
}

//Main function for the FEM powder layer consolidation
int main(int argc, char *argv[])
{
  printLine(__LINE__);
  
  //inititialize variables needed for calculations
  int myrank;
  int numranks;
  int i;
  ushort threadsCount = 0;
  int offset;
  //The only user input for the function is the desired number of threads
  threadsCount = atoi(argv[1]);

  printLine(__LINE__);
  
  //Commands to intialize MPI part
  MPI_Init(&argc, &argv);
  MPI_Status status; //Get status of MPI
  MPI_Request request; //Set up requests for MPI send and receive
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); //rank of the current process
  MPI_Comm_size(MPI_COMM_WORLD, &numranks); //Total Number of MPI ranks

  printLine(__LINE__);
  
  num_ElementsNodes(baseName, myrank);
  printLine(__LINE__);
  read_coordinates(baseName, myrank, nnodes);
  printLine(__LINE__);
  read_elements(baseName, myrank, nel);
  printLine(__LINE__);
  read_psi(baseName, myrank, nel);
  printLine(__LINE__);
  gol_runKernel(coordinates, nnodes, powderThick,
		Tol, elements,
		nel, &ID, threadsCount,
		&d, &a_bar);

  printLine(__LINE__);
  
  offset = offsetCalc(baseName, numranks, myrank);

  printLine(__LINE__);
  
  for(i=0;nel;i++)
  {
    elements[i] = elements[i] + offset[&myrank];
  }

  printLine(__LINE__);
  
  double timeStart; // start the timer for recording the processing time for the file writing
  if (myrank == 0)
    {
      printf("StartTimer\n");
      timeStart = MPI_Wtime();
    }

  //Write data to single file for use in post-processing of mesh consolidation
  MPI_File fh;
  int bufsize, nintsC, nintsE;
  //int buf[40000000];
  int fileLength = (sizeof(coordinates) + sizeof(elements))*numranks;
  bufsize = fileLength/numranks;
  nintsC = sizeof(coordinates)/sizeof(coordinates[0]);
  nintsE = sizeof(elements)/sizeof(elements[0]);

  printLine(__LINE__);
  
  MPI_File_open(MPI_COMM_WORLD,"Layer_270_09_10.vtk",MPI_MODE_RDWR |
		MPI_MODE_CREATE,
		MPI_INFO_NULL,&fh);
  MPI_File_write_at(fh,myrank*bufsize,coordinates,nintsC,MPI_FLOAT,&status);
  MPI_File_write_at(fh,myrank*bufsize+sizeof(coordinates),elements,nintsE,MPI_FLOAT,&status);
  MPI_File_close(&fh);

  printLine(__LINE__);
  
  double timeEnd;
  double time;
  //stop the timer and calculate the total run time for the file writing
  if (myrank == 0)
    {
      timeEnd = MPI_Wtime();
      time = timeEnd - timeStart;
      printf("Run time is %lf\n", time);
    }

  printLine(__LINE__);
  
  float coordBuf[sizeof(coordinates)];
  float elemBuf[sizeof(elements)];
  MPI_File_open(MPI_COMM_WORLD,"Layer_270_09_10.vtk",MPI_MODE_RDWR,
		MPI_INFO_NULL,&fh);
  MPI_File_read_at(fh,myrank*bufsize,coordBuf,nintsC,MPI_FLOAT,&status);
  MPI_File_read_at(fh,myrank*bufsize+sizeof(coordinates),elemBuf,nintsE,MPI_FLOAT,&status);
  MPI_File_close(&fh);

  printLine(__LINE__);
  
  MPI_Finalize();

  return false;
  //if(myrank == 0)
  //  {
  //    MPI_File_open(MPI_COMM_WORLD,"Layer_270_09_10.vtk",MPI_MODE_RDWR |
  //                MPI_MODE_CREATE | MPI_MODE_DELETE_ON_CLOSE,
  //		    MPI_INFO_NULL,&fh);
  //     char string = "# vtk DataFile Version 3.8\nVelocities, pressures and level set\nASCII\nDATASET UNSTRUCTURED_GRID\n";
  //  }
  
}

typedef unsigned long long ticks;
static __inline__ ticks getticks(void)
{
  unsigned int tbl, tbu0, tbu1;
  do {
      __asm__ __volatile__ ("mftbu %0" : "=r"(tbu0));
    __asm__ __volatile__ ("mftb %0" : "=r"(tbl));
    __asm__ __volatile__ ("mftbu %0" : "=r"(tbu1));
  } while (tbu0 != tbu1);
  return ((((unsigned long long)tbu0) << 32) | tbl);
}
