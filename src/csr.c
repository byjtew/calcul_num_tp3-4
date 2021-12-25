#include <float.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NULL_VALUE 0.0
#define TEST_DIMENSION 5000

//
double randrealf() {
    return ((double)rand() / RAND_MAX);
}

/**
 * @brief Compressed Sparse Row (CSR) Matrix type implementation
 */
typedef struct {
    double* AA;
    size_t AA_l;
    size_t IA_l;
    size_t *JA, *IA;
} mat_csr_t;

size_t get_mat_csr_nb_row(mat_csr_t* mat) { return mat->IA_l - 1; }
size_t get_mat_csr_nb_col(mat_csr_t* mat) {
    size_t maxi = 0;
    for (size_t i = 0; i < mat->AA_l; i++)
        if (maxi < mat->JA[i]) maxi = mat->JA[i];
    return maxi + 1;
}

void print_mat_csr_t(mat_csr_t* mat) {
    printf("\ndim: {%ld , %ld}\n", get_mat_csr_nb_row(mat),
           get_mat_csr_nb_col(mat));
    printf("AA: [");
    for (size_t i = 0; i < mat->AA_l; i++)
        printf(i == mat->AA_l - 1 ? "%.2lf" : "%.2lf, ", mat->AA[i]);
    printf("]\n");

    printf("JA: [");
    for (size_t i = 0; i < mat->AA_l; i++)
        printf(i == mat->AA_l - 1 ? "%ld" : "%ld, ", mat->JA[i]);
    printf("]\n");

    printf("IA: [");
    for (size_t i = 0; i < mat->IA_l; i++)
        printf(i == mat->IA_l - 1 ? "%ld" : "%ld, ", mat->IA[i]);
    printf("]\n\n");
}

double get_mat_element_at(mat_csr_t* mat, size_t row, size_t col) {
    // Bounds checking
    size_t nb_row = get_mat_csr_nb_row(mat);
    if (row >= nb_row) return NULL_VALUE;
    size_t nb_col = get_mat_csr_nb_col(mat);
    if (col >= nb_col) return NULL_VALUE;

    // No values for that line
    if (mat->IA[row + 1] - mat->IA[row] == 0) return NULL_VALUE;

    // Searching
    for (size_t i = mat->IA[row]; i < mat->IA[row + 1]; i++)
        if (mat->JA[i] == col) return mat->AA[i];

    return NULL_VALUE;
}

void print_mat_csr_t_full(mat_csr_t* mat) {
    printf("[\n");
    for (size_t l = 0; l < get_mat_csr_nb_row(mat); l++) {
        for (size_t c = 0; c < get_mat_csr_nb_col(mat); c++)
            printf("%.2lf, ", get_mat_element_at(mat, l, c));
        printf("\n");
    }
    printf("]\n");
}

void free_mat_csr_t(mat_csr_t* mat) {
    free(mat->AA);
    free(mat->IA);
    free(mat->JA);
}

mat_csr_t create_empty_mat_csr_t(size_t row) {
    mat_csr_t mat;
    mat.IA_l = row + 1;
    mat.IA = (size_t*)calloc(mat.IA_l, sizeof(size_t));
    mat.AA_l = 0;
    mat.JA = (size_t*)calloc(mat.AA_l, sizeof(size_t));
    mat.AA = (double*)calloc(mat.AA_l, sizeof(double));
    return mat;
}

mat_csr_t mat_csr_t_from_dense(double* mat_full, size_t row, size_t col) {
    mat_csr_t mat = create_empty_mat_csr_t(row);

    for (size_t r = 0; r < row; r++) {
        for (size_t c = 0; c < col; c++) {
            size_t element = mat_full[r * col + c];
            if (element != 0.0) {
                mat.AA_l += 1;
                mat.AA = (double*)realloc(mat.AA, mat.AA_l * sizeof(double));
                mat.AA[mat.AA_l - 1] = element;
                mat.JA = (size_t*)realloc(mat.JA, mat.AA_l * sizeof(size_t));
                mat.JA[mat.AA_l - 1] = c;
            }
        }
        mat.IA[r + 1] = mat.AA_l;
    }
    return mat;
}

/**
 * @brief Create a random mat_csr_t with a specified proportion of zeros
 *
 * @param n Dimension of the matrix
 * @param percent_of_zeros Obvious
 * @return mat_csr_t
 */
mat_csr_t create_random_csr_t(size_t n, double percent_of_zeros) {
    size_t nb_of_zeros = (size_t)(percent_of_zeros * n * n);
    double* full_mat = (double*)calloc(n * n, sizeof(double));
    for (size_t i = 0; i < nb_of_zeros; i++)
        full_mat[rand() % (n * n)] = randrealf();
    mat_csr_t mat = mat_csr_t_from_dense(full_mat, n, n);
    free(full_mat);
    return mat;
}

void print_vec(double* v, size_t n) {
    printf("[");
    for (size_t i = 0; i < n; i++)
        printf(i == n - 1 ? "%.2lf" : "%.2lf, ", v[i]);
    printf("]\n");
}

void mat_csr_mult_vec(mat_csr_t* mat, size_t n, double* vec, double* res_vec) {
    for (size_t r = 0; r < n; r++) {
        res_vec[r] = .0;
        for (size_t i = mat->IA[r]; i < mat->IA[r + 1]; i++)
            res_vec[r] += mat->AA[i] * vec[mat->JA[i]];
    }
}

int main(int argc, char** argv) {
    double start = 0, end = 0;
    double res_vector[TEST_DIMENSION] = {.0};
    double* res_vector_ptr = res_vector;

    // Saad CSR PART
    double saad_example[25] = {1.0, 0, 0, 2, 0,  3,  4, 0, 5, 0, 6, 0, 7,
                               8, 9, 0, 0, 10, 11, 0, 0, 0, 0, 0, 12};
    double* saad_ptr = saad_example;
    mat_csr_t saad_csr = mat_csr_t_from_dense(saad_ptr, 5, 5);
    printf("--- Saad example matrix ---\n");
    print_mat_csr_t(&saad_csr);
    print_mat_csr_t_full(&saad_csr);

    // PRODUCT CSR . VEC PART
    double two_vec[5] = {2.0, 2.0, 2.0, 2.0, 2.0};
    double* two_vec_ptr = two_vec;
    double res_vec[5] = {0.0};
    double* res_vec_ptr = res_vec;

    mat_csr_mult_vec(&saad_csr, 5, two_vec_ptr, res_vec_ptr);
    printf(
        "\n\n--- Multiplication of Saad example matrix by (2, ..., 2) "
        "vector:\n");
    print_vec(res_vec_ptr, 5);

    // RANDOM MATRIX PART
    printf("\n\n--- Random generated [%dx%d] matrix ---\n", TEST_DIMENSION, TEST_DIMENSION);
    srand(time(NULL));
    mat_csr_t random_csr = create_random_csr_t(TEST_DIMENSION, .75);

    // Question 1:
    printf("\n--- Question 1: Vector.ONE ---\n");
    double one_vec[TEST_DIMENSION];
    double* one_vec_ptr = one_vec;
    for (unsigned i = 0; i < TEST_DIMENSION; i++)
        one_vec[i] = i % 3 == 0 ? 1 : 0;

    // print_vec(one_vec_ptr, TEST_DIMENSION);
    for (size_t i = 0; i < TEST_DIMENSION; i++) res_vector_ptr[i] = .0;
    printf("\n-- Random [%dx%d] CSR matrix MULT Vector (1, 1, ...., 1)[%d] --\n",
           TEST_DIMENSION, TEST_DIMENSION, TEST_DIMENSION);
    start = omp_get_wtime();
    mat_csr_mult_vec(&random_csr, TEST_DIMENSION, one_vec_ptr, res_vector_ptr);
    // print_vec(res_vector_ptr, TEST_DIMENSION);
    end = omp_get_wtime();
    printf("-- Took %lf ms --\n", 1e3 * (end - start));

    // Question 2:
    printf("\n--- Question 2: Vector.ONE_ZERO_ZERO ---\n");
    double one_zero_zero_vec[TEST_DIMENSION];
    for (unsigned i = 0; i < TEST_DIMENSION; i++)
        one_zero_zero_vec[i] = i % 3 == 0 ? 1 : 0;
    double* one_zero_zero_vec_ptr = one_zero_zero_vec;
    // print_vec(one_zero_zero_vec_ptr, TEST_DIMENSION);
    for (size_t i = 0; i < TEST_DIMENSION; i++) res_vector_ptr[i] = .0;
    printf(
        "\n-- Random [%dx%d] CSR matrix MULT Vector (1, 0, 0, 1, 0, 0, ...., 1, "
        "0, 0)[%d] --\n",
        TEST_DIMENSION, TEST_DIMENSION, TEST_DIMENSION);
    start = omp_get_wtime();
    mat_csr_mult_vec(&random_csr, TEST_DIMENSION, one_zero_zero_vec_ptr,
                      res_vector_ptr);
    // print_vec(res_vector_ptr, TEST_DIMENSION);
    end = omp_get_wtime();
    printf("-- Took %lf ms --\n", 1e3 * (end - start));

    // End of the program
    free_mat_csr_t(&random_csr);
    free_mat_csr_t(&saad_csr);
}