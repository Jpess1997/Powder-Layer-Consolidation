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

unsigned char *coordinates;

unsigned char *elements;

size_t nel, nnodes, np, nzmax;

char basename;

unsigned char *psi;

unsigned char *LM;

unsigned char *irow;

unsigned char *icol;

float porosity = 0.8;

float powderThick = 28.5;

float Tol = 2.5;

basename = "/scratch/Powder-Layer-Consolidation/Layer_270_09_10/0/";

np = 4;

//extern functions from the gol-cuda file
extern void num_ElementsNodes(char basename, int myrank);

extern int offsetCalc(char basename, int numranks);

extern void read_coordinates(char basename, int myrank, size_t nnodes);

extern void read_elements(char basename, int myrank, size_t nel);

extern void read_psi(char basename, int myrank, size_t nel);

extern bool gol_runKernel(unsigned char coordinates, size_t nnodes, float powder_thick, float Tol, unsigned char elements,
			  size_t nel, ushort threadsCount);

//Main function for the FEM powder layer consolidation
int main(int argc, char *argv[])
{
  //inititialize variables needed for calculations
  int myrank;
  int numranks;
  int i;
  ushort threadsCount = 0;
  int offset;
  //The only user input for the function is the desired number of threads
  threadsCount = atoi(argv[1]);

  //Commands to intialize MPI part
  MPI_Init(&argc, &argv);
  MPI_Status status; //Get status of MPI
  MPI_Request request; //Set up requests for MPI send and receive
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); //rank of the current process
  MPI_Comm_size(MPI_COMM_WORLD, &numranks); //Total Number of MPI ranks

  num_ElementsNodes(basename, myrank);
  read_coordinates(basename, myrank, nnodes);
  read_elements(basename, myrank, nel);
  read_psi(basename, myrank, nel);
  gol_runKernel(coordinates, nnodes, powder_thick, Tol, elements,
		nel, ID, threadsCount);

  offset = offsetCalc(basename, numranks);

  for(i=0;nel;i++)
    {
      elements[i] = elements[i] + offset[myrank];
    }
