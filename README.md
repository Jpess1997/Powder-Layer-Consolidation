# Powder-Layer-Consolidation
Final project for parallel programming class Spring 2020. Adapts current code used to consolidate the powder layer as a post process for the thermal FEM simulation code.

This code is a prototype consolidation type that focuses on the implementation of parallel I/O, MPI, and CUDA for the use in a powder layer consolidation code for FEM simulation of an additive manufacturing process. The inputs for the code in the slurm file are first the number of threads per block for each rank, followed by the desired number of nnodes (the world size of each rank) in each rank. An example run would be:

gol-cuda-mpi-exe 1024 32768

This would be 1024 with 32768 nnodes.
