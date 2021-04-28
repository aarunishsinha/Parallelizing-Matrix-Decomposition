#include "mpi.h" 
#include <stdio.h> 
#include <string.h>
#include <stdlib.h>

void write_output(char fname[], double** arr, int n ){
	FILE *f = fopen(fname, "w");
	for( int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			fprintf(f, "%0.12f ", arr[i][j]);
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

int main(int argc, char **argv) 
{ 	
	char* end_pointer;
	int n = strtol(argv[1], &end_pointer, 10);
	char* filename = argv[2];
	
	FILE* f;
	double **inp = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0;i < n;i++){
        inp[i] = (double*)malloc(sizeof(double) * n);
	}

    if((f = fopen(filename, "r")) == NULL)
  	    exit(1);

    for(int jj = 0; jj < n; jj++){
        for(int ii = 0; ii<n; ii++){
            fscanf(f, "%lf", &inp[jj][ii]);
		}
	}
    fclose(f);
    const double **A = (const double **) inp;

    double **L = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0; i < n; i++){
        L[i] = (double*)malloc(sizeof(double) * n);
	}
	
	double **U = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0; i < n; i++){
        U[i] = (double*)malloc(sizeof(double)*n);
	}

    int size, my_rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    int i, j, k;
	double sum = 0;

	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			if(i % size == my_rank){
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L[i][k] * U[k][j];
				}
				L[i][j] = A[i][j] - sum;
			}
		}
		for(i = j; i < n; i++){
			MPI_Bcast(&(L[i][j]), 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
		}
		
		for (i = j; i < n; i++) {
			if(i % size == my_rank){
				sum = 0;
				for(k = 0; k < j; k++) {
					sum = sum + L[j][k] * U[k][i];
				}
				if (L[j][j] == 0) {
					exit(0);
				}
				U[j][i] = (A[j][i] - sum) / L[j][j];
			}
		}
		for(i = j; i < n; i++){
			MPI_Bcast(&(U[j][i]), 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
		}
		
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if(my_rank == 0){
		// printf("A\n");
		// print_output(A, n);
		
		char* ext = ".txt";
		char l_out_fname[50];
		strcpy(l_out_fname,"output_L_4_");
		const char * th = (const char *) argv[3];
		strncat(l_out_fname,th,strlen(th));
		strncat(l_out_fname,ext,4);

		write_output(l_out_fname,L,n);
		char u_out_fname[50];
		strcpy(u_out_fname,"output_U_4_");
		strncat(u_out_fname,th,strlen(th));
		strncat(u_out_fname,ext,4);
		write_output(u_out_fname,U,n);
	}


    MPI_Finalize();
    return 0;
 }