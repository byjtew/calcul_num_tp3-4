//#include "lib_poisson1D.h"
#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//
double randreal() {
    int a = rand() % 1000000;
    return ((double)a / 100000);
}

void print_tri_diag_matrix(double *A, unsigned n) {
    printf("[\n");
    for (unsigned i = 0; i < n; i++) printf("%.2lf, ", A[i]);
    printf("\n");
    for (unsigned i = n; i < 2 * n; i++) printf("%.2lf, ", A[i]);
    printf("\n");
    for (unsigned i = 2 * n; i < n * 3; i++) printf("%.2lf, ", A[i]);
    printf("]\n");
}

void tri_diag_zeros(double *A, unsigned n) {
    for (unsigned i = 0; i < n * 3; i++) A[i] = .0;
}

void tri_diag_randomised(double *A, unsigned n) {
    A[0] = .0;
    for (unsigned i = 1; i < n * 3 - 1; i++) A[i] = randreal();
    A[n * n - 1] = .0;
}


/**
 * @brief LU function for a square matrix
 *
 * @param A tri-diagonal matrix pointer
 * @param n leading dimension of matrix A
 */
void tri_diag_lu(double *A, unsigned n) {
    // Pointers to band index
    double *u = A, *d = A + n, *l = A + 2 * n;
    for (unsigned i = 1; i < n; i++) {
        printf("l[%d] /= d[%d]\n", i - 1, i - 1);
        *(l + i - 1) /= *(d + i - 1);
        printf("d[%d] -= ( l[%d] / d[%d]) * u[%d]\n", i, i - 1, i - 1, i - 1);
        *(d + i) -= *(l + i - 1) * (*(u + i));
    }
}

int main(int argc, char **argv) {
    unsigned dimension = 3;
    double *mat_A = calloc(dimension * 3, sizeof(double));

    // tri_diag_randomised(mat_A, dimension);
    mat_A[1] = -1;
    mat_A[2] = -6;
    mat_A[3] = 3;
    mat_A[4] = -2;
    mat_A[5] = 3;
    mat_A[6] = 2;
    mat_A[7] = 5;
    print_tri_diag_matrix(mat_A, dimension);

    tri_diag_lu(mat_A, dimension);
    print_tri_diag_matrix(mat_A, dimension);

    free(mat_A);
    return 0;
}