#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*--------------------------------------------------------
 * Procedure to allocate an n X n matrix and return the
 * pointer to the matrix.
 */
double** alloc_mtx(int n, int m)
{
	double** B;
	int     i;

	B = (double**)malloc(sizeof(double*) * n);
	for (i = 0; i < n; i++)
		B[i] = (double*)malloc(sizeof(double) * m);
	return (B);
}

/*-----------------------------------------------------
 * Procedure to allocate space for an n-dimensional
 * vector. Return the pointer to the vector.
 */
double* alloc_vec(int n)
{
	double* t;

	t = (double*)malloc(sizeof(double) * n);
	return(t);
}


/*--------------------------------------------------
 * Compute inner product <a, b>, where a and b are
 * n-dimensional vectors.
 */
double  inner_product(double* a, double* b, int n)
{
	int   i;
	double sum;

	sum = 0.0;
	for (i = 0; i < n; i++)
		sum += a[i] * b[i];
	return (sum);
}

/*---------------------------------------------------------
 * Compute a = A*b, where a and b are vectors and A an n x m
 * matrix.
 */
void mtx_vec_mult(double* a, double** A, double* b, int n, int m)
{
	int  i, j;

	for (i = 0; i < n; i++) {
		a[i] = 0.0;
		for (j = 0; j < m; j++)
			a[i] += A[i][j] * b[j];
	}
}

/*-------------------------------------------------------------
 * Normalize a vector(2norm). If the vector norm is 0, do nothing.
 */
void normalize_vec(double* a, int n)
{
	int    i;
	double sum;

	sum = 0.0;
	for (i = 0; i < n; i++)
		sum += a[i] * a[i];
	if (sum == 0.0) return;
	sum = sqrt(sum);
	for (i = 0; i < n; i++)
		a[i] = a[i] / sum;
}

/*--------------------------------------------------------
 * Procedure to compute the 2-norm of a vector.
 */
double vec_norm(double* a, int n)
{
	double   sum;
	int      i;

	sum = 0.0;
	for (i = 0; i < n; i++)
		sum += a[i] * a[i];
	sum = sqrt(sum);
	return(sum);
}

/*--------------------------------------------------------
 * Procedure to compute the 1-norm of a vector.
 */
double vec_1norm(double* a, int n)
{
	double   sum;
	int      i;

	sum = 0.0;
	for (i = 0; i < n; i++)
		sum += fabs(a[i]);
	return(sum);
}

/*------------------------------------------------
 * transposing A to A_T
 */
void transposing(double **A, double **A_T, int n, int m) {
	int i, j;
	
	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++) {
			A[i][j] = A_T[j][i];
		}
	}
}
/*-----------------------------------------------------
 * Copy a vector from source to destination.
 *    *src: source
 *    *dst: destination.
 */
void copy_vec(double* dst, double* src, int n)
{
	int   i;

	for (i = 0; i < n; i++) dst[i] = src[i];
}

/*-------------------------------------------------------------
 * Procedure to compute residual vector  r[] = x[] - u*y[].
 *    u: the eigen value.
 */
void comp_residual(double* r, double* x, double u, double* y, int n)
{
	int   i;

	for (i = 0; i < n; i++)
		r[i] = x[i] - u * y[i];
}

/*------------------------------------------------------
 * Procedure to print out a matrix.
 */
void print_mtx(double** A, int n, int m)
{
	int i, j;

	for (i = 0; i < n; i++) {
		fprintf(stderr, "\n");
		for (j = 0; j < m; j++)
			fprintf(stderr, "%lf ", A[i][j]);
	}
	fprintf(stderr, "\n--------------------------\n");
}


/*--------------------------------------------------------
 * Procedure to print out a vector.
 */
void print_vec(double* x, int n)
{
	int i;

	for (i = 0; i < n; i++) {
		fprintf(stderr, " %.10lf", x[i]);
	}
	fprintf(stderr, "\n--------------------------\n");
}