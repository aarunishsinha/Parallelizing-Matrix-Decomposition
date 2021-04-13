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
	int n = atoi(argv[1]);
	char* filename = argv[2];
	int num_threads = atoi(argv[3]);
	int strt = atoi(argv[4]);
	FILE* f;
	double inp[n][n];
	printf("Yo\n");
	if((f = fopen(filename, "r")) == NULL)
  	exit(1);
	printf("Yo1\n");
  // if(fscanf(f, "%f%f", &n, &n) != 2)
	// 	exit(1);
	// printf("Yo2\n");
  // if (height < 1 || height > MHEIGHT || width < 1 || width > MWIDTH)
  // 	exit(1);

  for(int jj=0; jj<n; jj++){
    for(int ii=0; ii<n; ii++){
      fscanf(f, "%lf", &inp[jj][ii]);
		}
	}
  fclose(f);
	printf("Y3\n");
  for(int jj=0; jj<n; jj++){
    for(int ii=0; ii<n; ii++)
      printf ("%0.12f", inp[jj][ii]);
    printf("\n");
  }
  return 0;
}
