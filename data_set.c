#define _CRT_SECURE_NO_WARNINGS
#include "function.h"

/*-------------------------
 * create 11 sample points with Horner' s algorithm.
 * y = p(x)= 1 + x + x^2 + x^3 + x^4 + x^5 + x^6 + x^7
 */
void create_sample_points(x_y_pair* sample_points, int n) {
    int i, j;
    double x;
    FILE* fp;

    fp = fopen("data_set12.txt", "w");

    for (j = 0, x = 2.0; j < n; j += 1, x += 0.2) {
        sample_points[j].x = x;
        sample_points[j].y = 1.0;
        for (i = 0; i < 7; i++) {
            sample_points[j].y = sample_points[j].y * x + 1;
        }
        fprintf(fp, "%lf %lf\n", sample_points[j].x, sample_points[j].y);
    }

    fclose(fp);
}


