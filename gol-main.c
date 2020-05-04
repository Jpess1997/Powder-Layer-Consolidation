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

char baseName[80]  = "Layer_270_09_10/0/";

float *psi;

float *LM;

int *irow;

int *icol;

float porosity = 0.8;

float powderThick = 28.5;

float Tol = 2.5;

size_t np = 4;

//extern functions from the gol-cuda file
extern void num_ElementsNodes(char baseName[80], int myrank, int nnodes);

//extern int offsetCalc(char baseName[80], int numranks, int myrank);

extern void read_coordinates(char baseName[80], int myrank, int nnodes);

extern void read_elements(char baseName[80], int myrank, int nel);

extern void read_psi(char baseName[80], int myrank, int nel);

extern bool gol_runKernel(float *coordinates, int nnodes, float powderThick,
			  float Tol, float *elements,
                          int nel, int **d_ID, ushort threadsCount,
			  float **d_d, float **d_a_bar);

extern void gol_freeData();

//This is an debugging function that prints out the most recent line the code gets to when running before an error occurs
void printLine(int line)
{
  char fileName[30]; //initializes file buffer
  snprintf(fileName,100,"errorAtLine.txt"); //created file name to be sent to
  FILE *fp;
  fp = fopen(fileName,"w+"); //opens file
  fprintf(fp,"Line is %d.\n",line); //prints the current line number to the file
  fclose(fp); //closes file
}

//function for the calculating the overhead for the code based on number of ticks
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

//Main function for the FEM powder layer consolidation
int main(int argc, char *argv[])
{
  printLine(__LINE__);

  //=============================================
  unsigned long long start = 0;
  unsigned long long finish = 0;

  start = getticks();
  //==============================================
  
  //inititialize variables needed for calculations
  int myrank;
  int numranks;
  int i;
  ushort threadsCount = 0;
  int offset;
  //The only user input for the function is the desired number of threads
  threadsCount = atoi(argv[1]);
  nnodes = atoi(argv[2]);

  printLine(__LINE__);
  
  //Commands to intialize MPI part
  MPI_Init(&argc, &argv);
  MPI_Status status; //Get status of MPI
  MPI_Request request; //Set up requests for MPI send and receive
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); //rank of the current process
  MPI_Comm_size(MPI_COMM_WORLD, &numranks); //Total Number of MPI ranks

  double timeStart; // start the timer for recording the time for CUDA and MPI calculations
  if (myrank == 0)
    {
      printf("StartTimer\n");
      timeStart = MPI_Wtime();
    }
  
  printLine(__LINE__);

  //run the code using CUDA and MPI for each rank
  num_ElementsNodes(baseName, myrank, nnodes);
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

  double timeEnd;
  double time;
  //stop the timer and calculate the total run time for the CUDA and MPI calculations
  if (myrank == 0)
    {
      timeEnd = MPI_Wtime();
      time = timeEnd - timeStart;
      printf("Run time for the cuda functions and MPI calculation is %lf\n", time);
    }
  
  printLine(__LINE__);

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
  nintsC = sizeof(coordinates)/sizeof(coordinates[0]); //set buffer for the coordinates to be written
  nintsE = sizeof(elements)/sizeof(elements[0]); //set buffer for the elements to be written

  printLine(__LINE__);
  
  MPI_File_open(MPI_COMM_WORLD,"Layer_270_09_10.vtk",MPI_MODE_RDWR |
		MPI_MODE_CREATE,
		MPI_INFO_NULL,&fh);
  MPI_File_write_at(fh,myrank*bufsize,coordinates,nintsC,MPI_FLOAT,&status);
  MPI_File_write_at(fh,myrank*bufsize+sizeof(coordinates),elements,nintsE,MPI_FLOAT,&status);
  MPI_File_close(&fh);

  printLine(__LINE__);
  
  //stop the timer and calculate the total run time for the file writing
  if (myrank == 0)
    {
      timeEnd = MPI_Wtime();
      time = timeEnd - timeStart;
      printf("Run time for the MPI I/O file writing is %lf\n", time);
    }

  printLine(__LINE__);

  //start timer for calcualting time to read in the data
  if (myrank == 0)
    {
      printf("StartTimer\n");
      timeStart = MPI_Wtime();
    }
  
  float coordBuf[sizeof(coordinates)]; //buffer for reading coordinates
  float elemBuf[sizeof(elements)]; //buffer for reading elements
  MPI_File_open(MPI_COMM_WORLD,"Layer_270_09_10.vtk",MPI_MODE_RDWR,
		MPI_INFO_NULL,&fh);
  MPI_File_read_at(fh,myrank*bufsize,coordBuf,nintsC,MPI_FLOAT,&status);
  MPI_File_read_at(fh,myrank*bufsize+sizeof(coordinates),elemBuf,nintsE,MPI_FLOAT,&status);
  MPI_File_close(&fh);

  printLine(__LINE__);

   //stop the timer and calculate the total run time for the file reading
  if (myrank == 0)
    {
      timeEnd = MPI_Wtime();
      time = timeEnd - timeStart;
      printf("Run time for the MPI I/O file reading is %lf\n", time);
    }
  
  MPI_Finalize();

  //============================================================================================
  //usleep(10000000); // 10 seconds sleep

  finish = getticks();

  printf("10 second usleep: finish(%llu) - start(%llu) = %llu \n", finish, start, (finish-start));
  //============================================================================================
  
  return false;

}
