#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
double **L_master;
double **U_master;
void write_complete(double **L,double **U,int n){
	for(int i =0;i<n;i++){
		for(int j=0;j<n;j++){
			L_master[i][j]=L[i][j];
			U_master[i][j]=U[i][j];
		}
	}
}
void write_L(double **L,int n){
	for(int i=0;i<n;i++){
		L_master[i][i]=L[i][i];
	}
}
void write_U1(int n){
	for(int i=0;i<n;i++){
		U_master[i][i]=1;
	}
}
void write_U(double **U,int n){
	for(int j=0;j<n;j++){
		for(int i = j; i < n; i++){
			U_master[j][i]=U[j][i];
		}
	}
}
void crout(double const **A,double **L,double **U, int n, int argc, char* argv[]){
	int i, j, k, comm_sz, my_rank;
	double sum = 0;
	double start=0,end=0;
	double **L_dash = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0; i < n; i++){
        L_dash[i] = (double*)malloc(sizeof(double) * n);
	}
	double **U_dash = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0; i < n; i++){
        U_dash[i] = (double*)malloc(sizeof(double)*n);
	}
	MPI_Init(&argc,&argv);
	// MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);
	start = MPI_Wtime();
	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}
	write_U1(n);
	if(my_rank==1){
		for(j=0;j<n;j++){
			sum = 0;
			for (k = 0; k < j; k++) {
					sum = sum + L[j][k] * U[k][j];
			}
			L[j][j] = A[j][j] - sum;
			for (i = j; i < n; i++) {
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
		for(j=0;j<n;j++){
			MPI_Send(&L[j][j],1,MPI_DOUBLE,0,j,MPI_COMM_WORLD);
		}
		for(j=0;j<n;j++){
			for(i = j; i < n; i++){
				MPI_Send(&U[j][i],1,MPI_DOUBLE,0,(n+i+n*j),MPI_COMM_WORLD);
			}
		}
	}
	else if(my_rank==0){
		for(j=0;j<n;j++){
			for (i = j + 1; i < n; i++) {
				sum = 0;
				for (k = 0; k < j; k++) {
						sum = sum + L[i][k] * U[k][j];
				}
				L[i][j] = A[i][j] - sum;
			}
		}
		for(j=0;j<n;j++){
			MPI_Recv(&L_dash[j][j],1,MPI_DOUBLE,1,j,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		}
		for(j=0;j<n;j++){
			for(i = j; i < n; i++){
				MPI_Recv(&U_dash[j][i],1,MPI_DOUBLE,1,(n+i+n*j),MPI_COMM_WORLD,MPI_STATUS_IGNORE);
			}
		}
		write_complete(L,U,n);
		write_L(L_dash,n);
		write_U(U_dash,n);
	}
	end = MPI_Wtime();
	MPI_Finalize();
	printf("Work took %f seconds\n", end - start);
}
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

int main(int argc, char* argv[]){
	// double start = 0, end = 0;
    char* end_pointer;
	int n = strtol(argv[1], &end_pointer, 10);
	char* filename = argv[2];
	int num_threads = strtol(argv[3], &end_pointer, 10);
	int strt = strtol(argv[4], &end_pointer, 10);
	FILE* f;
	double **inp = (double**)malloc(n * (sizeof(double*)));
	for(int i=0;i<n;i++){
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
	// Printing the matrix
    // for(int jj=0; jj<n; jj++){
    //   for(int ii=0; ii<n; ii++)
    //     printf ("%0.12f", A[jj][ii]);
    //   printf("\n");
    // }
	L_master = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0; i < n; i++){
		L_master[i] = (double*)malloc(sizeof(double) * n);
	}
	U_master = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0; i < n; i++){
		U_master[i] = (double*)malloc(sizeof(double)*n);
	}
	double **L = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0; i < n; i++){
        L[i] = (double*)malloc(sizeof(double) * n);
	}
	double **U = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0; i < n; i++){
        U[i] = (double*)malloc(sizeof(double)*n);
	}
	crout(A,L,U,n,argc,argv);
	char* ext = ".txt";
	char l_out_fname[50];
	strcpy(l_out_fname,"output_L_4_");
	const char * th = (const char *) argv[3];
	strncat(l_out_fname,th,strlen(th));
	strncat(l_out_fname,ext,4);

	write_output(l_out_fname,L_master,n);
	char u_out_fname[50];
	strcpy(u_out_fname,"output_U_4_");
	strncat(u_out_fname,th,strlen(th));
	strncat(u_out_fname,ext,4);
	write_output(u_out_fname,U_master,n);
	// printf("Work took %f seconds\n", end - start);
  return 0;
}
