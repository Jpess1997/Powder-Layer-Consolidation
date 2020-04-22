//Main function code to run on the CPU
//Jacob Pessin & Alex Vest, April 2020 Parallel Programming final project
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

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

basename = "/Users/Jacob/Desktop/First_Year_PhD/Research/Layer_270_09_10/0/";

np = 4;


