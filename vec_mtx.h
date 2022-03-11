/*--------------------------------------------------------
 * Procedure to allocate an n X n matrix and return the
 * pointer to the matrix.
 */
double** alloc_mtx(int n, int m);

/*-----------------------------------------------------
 * Procedure to allocate space for an n-dimensional
 * vector. Return the pointer to the vector.
 */
double* alloc_vec(int n);

/*--------------------------------------------------
 * Compute inner product <a, b>, where a and b are
 * n-dimensional vectors.
 */
double  inner_product(double* a, double* b, int n);

/*---------------------------------------------------------
 * Compute a = A*b, where a and b are vectors and A an n x m
 * matrix.
 */
void mtx_vec_mult(double* a, double** A, double* b, int n, int m);

/*------------------------------------------------
 * transposing A to A_T
 */
void transposing(double** A, double** A_T, int n, int m);

/*-------------------------------------------------------------
 * Normalize a vector. If the vector norm is 0, do nothing.
 */
void normalize_vec(double* a, int n);

/*--------------------------------------------------------
 * Procedure to compute the 2-norm of a vector.
 */
double vec_norm(double* a, int n);

/*--------------------------------------------------------
 * Procedure to compute the 1-norm of a vector.
 */
double vec_1norm(double* a, int n);

/*-----------------------------------------------------
 * Copy a vector from source to destination.
 *    *src: source
 *    *dst: destination.
 */
void copy_vec(double* src, double* dst, int n);

/*-------------------------------------------------------------
 * Procedure to compute residual vector  r[] = x[] - u*y[].
 *    u: the eigen value.
 */
void comp_residual(double* r, double* x, double u, double* y, int n);

/*------------------------------------------------------
 * Procedure to print out a matrix.
 */
void print_mtx(double** A, int n, int m);

/*--------------------------------------------------------
 * Procedure to print out a vector.
 */
void print_vec(double* x, int n); 
