#!/bin/bash
gcc -O0 -o omp_strats -fopenmp omp_strats.c
mpicc -o mpi_strat mpi_strat.c
