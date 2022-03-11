#define _CRT_SECURE_NO_WARNINGS
#include "function.h"
#include "vec_mtx.h"

/*--------------------------------------------------------
 * Initialize matrix and vector.
 *   A: the coef.mtx, 
 *   y: the right hand side.
 *   sample_points: sample_points,
 *   n: numbers of rows of matrix = 11
 *   m: numbers of columns of matrix = 8 
 */
void initial_mtx_n_vec(double** A, double** A_T, double* y, x_y_pair* sample_points, int n, int m) {
    double tmp;
    for (int i = 0; i < n; i++) {
        y[i] = sample_points[i].y;
        A[i][0] = 1.0;
        tmp = sample_points[i].x;
        for (int j = 1; j < m; j++) {
            A[i][j] = tmp;
            tmp *= sample_points[i].x;
        }
    }
    transposing(A_T, A, m, n);
}

/*------------------------------------------------------
 * Create a symmetric linear system.
 *   A: the coef. mtx,
 *   b: the right hand side.
 *   n: dimension of matrix.
 */
void symmetric_linear_system(double** A, double* b, int n)
{
    int  i, j;
    int  t;

    srand(0);
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            t = rand() % 10;
            A[i][j] = A[j][i] = (double)(t + 1.0);
        }
    }
    for (i = 0; i < n; i++) {
        b[i] = 0.0;
        for (j = 0; j < n; j++) b[i] += A[i][j];
    }
}

/*-------------------------------------------------------
 * B = A_T * A.
 */
void new_mtx_n_vec(double** B, double** A, double** A_T, double* y, double* d, int n, int m) {
    for (int i = 0; i < m; i++) {
        d[i] = 0.0;
        for (int j = 0; j < n; j++)
            d[i] += A_T[i][j] * y[j];
    }
    for (int x = 0; x < m; x++) {
        for (int y = 0; y < m; y++) {
            B[x][y] = 0.0;
            for (int i = 0; i < n; i++) {
                B[x][y] += A_T[x][i] * A[i][y];
            }
        }
    }
}


/**************************************************************
 * Perform forward elimination for a linear system Ax=b.
 * Partial pivoting is adopted.
 */
void gauss_elm(double** A, double* b, int n)
{
    int     p, i, j, k;
    double  maxEntry, t, r;

    for (i = 0; i < n - 1; i++) {
        // Partial pivoting
        maxEntry = fabs(A[i][i]);
        p = i;
        for (k = i; k < n; k++)
            if (fabs(A[k][i]) > maxEntry) {
                p = k;
                maxEntry = fabs(A[k][i]);
            }
        if (i != p) {
            t = b[p];
            b[p] = b[i];
            b[i] = t;
            for (j = i; j < n; j++) {
                t = A[i][j];
                A[i][j] = A[p][j];
                A[p][j] = t;
            }
        }
        //Forward elimination.
        for (k = i + 1; k < n; k++) {
            if (A[k][i] == 0.0) continue;

            r = A[k][i] / A[i][i];
            A[k][i] = 0.0;
            for (j = i + 1; j < n; j++) {
                A[k][j] -= r * A[i][j];
                //printf_s("A[%d][%d]: %lf\t", k, j, A[k][j]);
            }
            b[k] = b[k] - r * b[i];
        }
    }
}

/*
 * free matrix
*/
void free_mtx(double** B, int n, int m)
{
    int i;
    for (i = n-1; i >= 0; i--)
        free(B[i]);
}

/*
 * compute residual of the least squares.
 * r[] = b[] - A[][]*x[]
               (   t[]  )
 */
void residual(double* r, double* b, double* t, int n) {
    for (int i = 0; i < n; i++) {
        r[i] = b[i] - t[i];
    }
}

/*
 * printf points(x,y) to file.
 */
void printf_points_to_file(int num, x_y_pair* sample_points, double* t, int n) {
    
    FILE* fp;
    int i = 0;
    if(num == 1)
        fp = fopen("result1_points.txt", "w");
    else if(num == 2)
        fp = fopen("result2_points.txt", "w");
    else 
        fp = fopen("result3_points.txt", "w");

    for (i = 0; i < n; i++) {
        fprintf(fp, "%lf %lf\n", sample_points[i].x, t[i]);
    }
    fclose(fp);
}

int main() {

    x_y_pair sample_points[11];
    
    /*--- Declare the coef. matrix, the unkown vec. and the rhs.. ---*/
    double** A, **A_T, * c, * y, ** B, * d, * r, * t;
    double two_norm_residual, one_norm_residual, tmp;
    int n = 11, m = 8, i, j; // dimension of the system.
    FILE* fp;
    
    create_sample_points(sample_points, n);

    /*
    A*c = y => A_T*A*c = A_T*y => B*c = d
    t = A*c;  r = y-t; 
    */
    A = alloc_mtx(n, m);
    A_T = alloc_mtx(m, n);
    c = alloc_vec(m);
    y = alloc_vec(n);
    B = alloc_mtx(m, m);
    d = alloc_vec(m);
    r = alloc_vec(n); 
    t = alloc_vec(n);


    /*
     * initialize matrix A and vector y
     * and make new system
    */
    initial_mtx_n_vec(A, A_T, y, sample_points, n, m);
    new_mtx_n_vec(B, A, A_T, y, d, n, m);

    /* ----------------------------------------------------
     * result1
     * gaussian elimination forwording
     * backwording
    */
    /*
    //print_mtx(B, m);
    gauss_elm(B, d, m);
    print_mtx(B, m, m);
    back_substitute(B, c, d, m);
    print_vec(c, m);
    mtx_vec_mult(t, A, c, n, m);
    residual(r, y, t, n, m);
    two_norm_residual = vec_norm(r, n);
    one_norm_residual = vec_1norm(r, n);
    printf_points_to_file(1, sample_points, t, n);
    printf("two_norm_residual: %lf\n", two_norm_residual);
    printf("one_norm_residual: %lf\n", one_norm_residual);
    */

    /*-------------------------------------------------------------- 
     * result2
     * QR_new_system
    */
    /*
    
    QR_solver(B, c, d, m, m);
    print_vec(c, m);
    mtx_vec_mult(t, A, c, n, m);
    //print_vec(t, n);
    residual(r, y, t, n, m);
    //print_vec(r, n);
    two_norm_residual = vec_norm(r, n);
    one_norm_residual = vec_1norm(r, n);
    printf_points_to_file(2, sample_points, t, n);
    printf("two_norm_residual: %.10lf\n", two_norm_residual);
    printf("one_norm_residual: %.10lf\n", one_norm_residual);
    
    */

    /* --------------------------------------------------------
     * result3
     * QR_original_system
    */
    
    QR_solver(A, c, y, n, m);
    print_vec(c, m);

    //resuming A and y
    initial_mtx_n_vec(A, A_T, y, sample_points, n, m);
    //print_mtx(A, n, m);
    mtx_vec_mult(t, A, c, n, m);
    print_vec(t, n);
    residual(r, y, t, n);
    print_vec(r, n);
    two_norm_residual = vec_norm(r, n);
    one_norm_residual = vec_1norm(r, n);
    printf_points_to_file(3, sample_points, t, n);
    printf("two_norm_residual: %.10lf\n", two_norm_residual);
    printf("one_norm_redisual: %.10lf\n", one_norm_residual);
    

    

 
    free(c);
    free(y);
    free(d);
    free(r);
    free_mtx(B, m, m);
    free_mtx(A, n, m);
    free_mtx(A_T, m, n);

    return 0;
}