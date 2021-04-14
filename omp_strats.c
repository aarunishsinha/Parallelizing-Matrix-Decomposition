#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
void crout(double const **A, double **L, double **U, int n) {
	int i, j, k;
	double sum = 0;

	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}

	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			sum = 0;
			for (k = 0; k < j; k++) {
				sum = sum + L[i][k] * U[k][j];
			}
			L[i][j] = A[i][j] - sum;
		}

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
}
void s2_crout(double const **A, double **L, double **U, int n, int num_threads) {
	int i, j, k;
	double sum = 0;
	// omp_set_num_threads(num_threads);
	for (i = 0; i < n; i++) {
		U[i][i] = 1;
	}
	for (j = 0; j < n; j++) {
		#pragma omp parallel num_threads(num_threads) private(sum)
		{
			#pragma omp sections private(sum)
			{
				#pragma omp section
				{
					for (i = j+1; i < n; i++) {
						sum = 0;
							for (k = 0; k < j; k++) {
								sum = sum + L[i][k] * U[k][j];
							}
						L[i][j] = A[i][j] - sum;
					}
				}
			}
			#pragma omp sections private(sum)
			{
				#pragma omp section
				{
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
			}
		}
	}
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
	double start,end;
	int n = atoi(argv[1]);
	char* filename = argv[2];
	int num_threads = atoi(argv[3]);
	int strt = atoi(argv[4]);
	FILE* f;
	double **inp=(double**)malloc(n*(sizeof(double*)));
	for(int i=0;i<n;i++){
    inp[i]=(double*)malloc(sizeof(double)*n);
	}
	if((f = fopen(filename, "r")) == NULL)
  	exit(1);

  for(int jj=0; jj<n; jj++){
    for(int ii=0; ii<n; ii++){
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
	double **L=(double**)malloc(n*(sizeof(double*)));
	for(int i=0;i<n;i++){
    L[i]=(double*)malloc(sizeof(double)*n);
	}
	double **U=(double**)malloc(n*(sizeof(double*)));
	for(int i=0;i<n;i++){
    U[i]=(double*)malloc(sizeof(double)*n);
	}
	if(strt==0){
		start = omp_get_wtime();
		crout(A,L,U,n);
		end = omp_get_wtime();
		char* ext = ".txt";
		char l_out_fname[50];
		strcpy(l_out_fname,"output_L_0_");
		const char * th = (const char *) argv[3];
		strncat(l_out_fname,th,strlen(th));
		strncat(l_out_fname,ext,4);

		write_output(l_out_fname,L,n);
		char u_out_fname[50];
		strcpy(u_out_fname,"output_U_0_");
		strncat(u_out_fname,th,strlen(th));
		strncat(u_out_fname,ext,4);
		write_output(u_out_fname,U,n);
	}
	else if(strt==2){
		start = omp_get_wtime();
		s2_crout(A,L,U,n,num_threads);
		end = omp_get_wtime();
		char* ext = ".txt";
		char l_out_fname[50];
		strcpy(l_out_fname,"output_L_2_");
		const char * th = (const char *) argv[3];
		strncat(l_out_fname,th,strlen(th));
		strncat(l_out_fname,ext,4);

		write_output(l_out_fname,L,n);
		char u_out_fname[50];
		strcpy(u_out_fname,"output_U_2_");
		strncat(u_out_fname,th,strlen(th));
		strncat(u_out_fname,ext,4);
		write_output(u_out_fname,U,n);
	}

	printf("Work took %f seconds\n", end - start);
  return 0;
}
