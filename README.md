# Parallelizing-Matrix-Decomposition
COL380 - Parallel and Distributed Programming

## Crout Matrix Decomposition
### Theory
The Crout matrix decomposition is an LU decomposition that decomposes a matrix into a lower triangular matrix (L), an upper triangular matrix (U) and, although not always needed, a permutation matrix (P). It was developed by Prescott Durand Crout. Crout method returns a lower triangular matrix and a unit upper triangular matrix.\
So, if a matrix decomposition of a matrix A is such that:\
A = LDU\
being L a unit lower triangular matrix, D a diagonal matrix and U a unit upper triangular matrix, then Crout's method produces\
A = (LD)U = LU\
### Sequential Program
```C
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
```
Reference: [Wikipedia](https://en.wikipedia.org/wiki/Crout_matrix_decomposition)

## Methods
### Strategy 1
In this strategy we were supposed to use the ```parallel for``` construct of ```OpenMP```.
### Strategy 2
In this strategy we were supposed to use the ```parallel sections``` construct of ```OpenMP```.
### Strategy 3
In this strategy we were supposed to use both ```parallel for``` and ```parallel sections``` construct of ```OpenMP```.
### Strategy 4
In this strategy we were supposed to write an ```MPI``` version that solves the problem in a distributed manner.

**Contraint**: Do not use `reduction` or `atomic` clauses in any implementation.

For more details and analysis check [Report](https://github.com/aarunishsinha/Parallelizing-Matrix-Decomposition/blob/master/Report.pdf)
