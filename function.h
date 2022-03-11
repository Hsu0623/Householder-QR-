#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct x_y_pair x_y_pair;
struct x_y_pair {
    double x;
    double y;
};

/*-----------------------------------------------------------
 * create 11 sample points with Horner' s algorithm.
 * y = p(x)= 1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7 
 */
void create_sample_points(x_y_pair* sample_points, int n);


/*------------------------------------------------------------
 * Procedure to solve a lower triangular system
 *    L*x = b.
 */
void back_substitute(double** U, double* x, double* b, int n);


/**************************************************************
 * Procedure to solve a linear system by using QR-decompsotion.
 * Algm:
 *    1. convert Ax=b into an upper triangular system.
 *    2. Solve the upper trianglular system by using
 *       backward substitution.
 */
void QR_solver(double** A, double* x, double* b, int n, int m);

