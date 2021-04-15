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

void s1_crout(double const **A, double **L, double **U, int n, int num_threads) {
    int i, j, k;
    double sum = 0;

    #pragma omp parallel num_threads(num_threads) private(sum, i, j, k)
    {
        #pragma omp for
        for (i = 0; i < n; i++) {
            U[i][i] = 1;
        }

        for (j = 0; j < n; j++) {
            #pragma omp for
            for (i = j; i < n; i++) {
                sum = 0;
                for (k = 0; k < j; k++) {
                    sum = sum + L[i][k] * U[k][j];
                }
                L[i][j] = A[i][j] - sum;
            }
            #pragma omp for
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

void s2_crout(double const **A, double **L, double **U, int n, int num_threads) {
    int i, j, k;
    // omp_set_num_threads(num_threads);
    double sum = 0;
    for (i = 0; i < n; i++) {
        U[i][i] = 1;
    }
    for (j = 0; j < n; j++) {
        #pragma omp parallel sections num_threads(num_threads) private(sum, i, k)
        {
            #pragma omp section
            {
                for (i = j + 1; i < n; i++) {
                    sum = 0;
                    for (k = 0; k < j; k++) {
                        sum = sum + L[i][k] * U[k][j];
                    }
                    L[i][j] = A[i][j] - sum;
                }
            }
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

void s3_section_1(double const **A, double **L, double **U, int n, int num_threads, int j){
    #pragma omp parallel for num_threads(num_threads)
    for (int i = j + 1; i < n; i++) {
        double sum = 0;
        for (int k = 0; k < j; k++) {
            sum = sum + L[i][k] * U[k][j];
        }
        L[i][j] = A[i][j] - sum;
    }
}

void s3_section_2(double const **A, double **L, double **U, int n, int num_threads, int j){
    double sum = 0;
    for (int k = 0; k < j; k++) {
        sum = sum + L[j][k] * U[k][j];
    }
    L[j][j] = A[j][j] - sum;
    #pragma omp parallel for private(sum) num_threads(num_threads)
    for (int i = j; i < n; i++) {
        sum = 0;
        for(int k = 0; k < j; k++) {
            sum = sum + L[j][k] * U[k][i];
        }
        if (L[j][j] == 0) {
            exit(0);
        }
        U[j][i] = (A[j][i] - sum) / L[j][j];
    }
}

void s3_crout(double const **A, double **L, double **U, int n, int num_threads) {
    int i, j, k;
    // double sum = 0;
    omp_set_nested(1);

    #pragma omp parallel for private(i) num_threads(num_threads)
    for (i = 0; i < n; i++) {
        U[i][i] = 1;
    }

    for (j = 0; j < n; j++) {
        #pragma omp parallel sections num_threads(2)
        {
            #pragma omp section
            {
                s3_section_1(A, L, U, n, num_threads/2, j);
            }

            #pragma omp section
            {
                s3_section_2(A, L ,U, n, num_threads - num_threads/2, j);
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
	double start = 0, end = 0;
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
	double **L = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0; i < n; i++){
        L[i] = (double*)malloc(sizeof(double) * n);
	}
	double **U = (double**)malloc(n * (sizeof(double*)));
	for(int i = 0; i < n; i++){
        U[i] = (double*)malloc(sizeof(double)*n);
	}
    if(strt == 0){
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
    else if(strt == 1){
        start = omp_get_wtime();
        s1_crout(A,L,U,n,num_threads);
        end = omp_get_wtime();
        char* ext = ".txt";
        char l_out_fname[50];
        strcpy(l_out_fname,"output_L_1_");
        const char * th = (const char *) argv[3];
        strncat(l_out_fname,th,strlen(th));
        strncat(l_out_fname,ext,4);

        write_output(l_out_fname,L,n);
        char u_out_fname[50];
        strcpy(u_out_fname,"output_U_1_");
        strncat(u_out_fname,th,strlen(th));
        strncat(u_out_fname,ext,4);
        write_output(u_out_fname,U,n);
    }
    else if(strt == 2){
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
    else if(strt == 3){
        start = omp_get_wtime();
        s3_crout(A,L,U,n,num_threads);
        end = omp_get_wtime();
        char* ext = ".txt";
        char l_out_fname[50];
        strcpy(l_out_fname,"output_L_3_");
        const char * th = (const char *) argv[3];
        strncat(l_out_fname,th,strlen(th));
        strncat(l_out_fname,ext,4);

        write_output(l_out_fname,L,n);
        char u_out_fname[50];
        strcpy(u_out_fname,"output_U_3_");
        strncat(u_out_fname,th,strlen(th));
        strncat(u_out_fname,ext,4);
        write_output(u_out_fname,U,n);
    }
    else if(strt == 4){
        printf("Not Implemented\n");
        return 0;
    }
    else{
        exit(EXIT_FAILURE);
    }
	printf("Work took %f seconds\n", end - start);
  return 0;
}
