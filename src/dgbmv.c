#include <unistd.h>

#include "blaslapack_headers.h"
#include "math.h"

typedef double real;
const size_t real_size = sizeof(real);

real randreal() { return (real)rand() / (real)RAND_MAX; }

void print_tri_diag_matrix(real *A_tri_diag, size_t n) {
    printf("[\n");
    for (size_t i = 0; i < n; i++) printf("%lf ", A_tri_diag[i]);
    printf("\n");
    for (size_t i = n; i < 2 * n; i++) printf("%lf ", A_tri_diag[i]);
    printf("\n");
    for (size_t i = 2 * n; i < n * 3; i++) printf("%lf ", A_tri_diag[i]);
    printf("\n]\n");
}
void print_dense_matrix(real *A_dense, size_t n) {
    printf("[\n");
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) printf("%lf ", A_dense[i * n + j]);
        printf("\n");
    }
    printf("]\n");
}
void print_vector(real *vector, size_t n) {
    printf("[ ");
    for (size_t i = 0; i < n; i++) printf("%lf ", vector[i]);
    printf("]\n");
}

void tri_diag_matrix_from_dense_rowMajor(real *A_tri_diag, real *A_dense,
                                         size_t n) {
    real *ub_ptr = A_tri_diag + 1, *d_ptr = A_tri_diag + n,
         *lb_ptr = A_tri_diag + 2 * n;
    for (size_t i = 0; i < n - 1; i++) {
        *(ub_ptr + i) = A_dense[i * n + 1 + i];
        *(d_ptr + i) = A_dense[i * n + i];
        *(lb_ptr + i) = A_dense[n + i * n + i];
    }
    *(d_ptr + n - 1) = A_dense[n * n - 1];
    print_tri_diag_matrix(A_tri_diag, n);
}

void tri_diag_matrix_from_dense_colMajor(real *A_tri_diag, real *A_dense,
                                         size_t n) {
    for (unsigned i = 0; i < n - 1; i++) {
        A_tri_diag[(i + 1) * 3] = A_dense[i * n + (i + 1)];
        A_tri_diag[i * 3 + 1] = A_dense[i * n + i];
        A_tri_diag[i * 3 + 2] = A_dense[(i + 1) * n + i];
    }
    A_tri_diag[(n - 1) * 3 + 1] = A_dense[n * n - 1];
    print_tri_diag_matrix(A_tri_diag, n);
}

void tri_diag_matrix_from_dense(CBLAS_ORDER major, real *A_tri_diag,
                                real *A_dense, size_t n) {
    if (major == CblasColMajor)
        tri_diag_matrix_from_dense_colMajor(A_tri_diag, A_dense, n);
    else
        tri_diag_matrix_from_dense_rowMajor(A_tri_diag, A_dense, n);
}

void random_dense_matrix(size_t n, real *A_dense) {
    for (size_t i = 0; i < n - 1; i++) {
        A_dense[i * n + (i + 1)] = (real)rand() / (real)RAND_MAX;
        A_dense[i * n + i] = (real)rand() / (real)RAND_MAX;
        A_dense[(i + 1) * n + i] = (real)rand() / (real)RAND_MAX;
    }
    A_dense[n * n - 1] = randreal();

    print_dense_matrix(A_dense, n);
}

int main(int argc, char *argv[]) {
    int la = 6;
    int ku = 1, kl = 1, kv = 0, lab = 0;
    CBLAS_ORDER blas_order = CblasColMajor;

    lab = kv + kl + ku + 1;

    printf(
        "=== Generation of a random (%dx%d) matrix (dense & tri_diagonal "
        "formats) ===\n\n",
        la, la);

    real *A_tri_diag = calloc(la * lab, real_size);
    real *x_tri_diag = calloc(la, real_size);
    real *y_tri_diag = calloc(la, real_size);

    real *A_dense = calloc(la * la, real_size);
    real *x_dense = calloc(la, real_size);
    real *y_dense = calloc(la, real_size);

    random_dense_matrix(la, A_dense);
    tri_diag_matrix_from_dense(blas_order, A_tri_diag, A_dense, la);

    for (unsigned i = 0; i < la; i++) {
        real tmp_x = randreal();
        x_tri_diag[i] = tmp_x;
        x_dense[i] = tmp_x;

        real tmp_y = randreal();
        y_dense[i] = tmp_y;
        y_tri_diag[i] = tmp_y;
    }

    real alpha = randreal();
    real beta = randreal();

    // y_tri_diag := alpha * A_tri_diag * x_tri_diag + beta * y_tri_diag
    cblas_dgbmv(blas_order, CblasNoTrans, la, la, kl, ku, alpha, A_tri_diag,
                lab, x_tri_diag, 1, beta, y_tri_diag, 1);

    printf(
        "\n=== Check result by comparing DGEMV result and DGBMV result "
        "===\n\n");

    cblas_dgemv(blas_order, CblasNoTrans, la, la, alpha, A_dense, la, x_dense,
                1, beta, y_dense, 1);

    print_vector(y_tri_diag, la);
    print_vector(y_dense, la);

    real diff = .0;
    for (size_t i = 0; i < la; i++) diff += fabs(y_dense[i] - y_tri_diag[i]);
    diff /= (real)la;

    if (diff <= 1e-15)
        printf(
            "Les vecteurs y_tri_diag et y_dense sont egaux: mean absolute "
            "difference = %lf\n",
            diff);
    else
        printf(
            "Les vecteurs y_tri_diag et y_dense ne sont pas egaux: mean "
            "absolute difference = %lf\n",
            diff);

    real norme_y_tri_diag = cblas_dnrm2(la, y_tri_diag, 1);
    real norme_y_dense = cblas_dnrm2(la, y_dense, 1);

    printf("Norme y_tri_diag : %lf\n", norme_y_tri_diag);
    printf("Norme y_dense : %lf\n", norme_y_dense);
    fflush(stdout);

    printf("\n--------- End -----------\n");

    free(A_tri_diag);
    free(x_tri_diag);
    free(y_tri_diag);
    free(A_dense);
    free(x_dense);
    free(y_dense);

    return 0;
}
