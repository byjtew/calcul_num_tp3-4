#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NULL_VALUE 0.0

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
        printf(i == mat->AA_l - 1 ? "%lf" : "%lf, ", mat->AA[i]);
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
            printf("%lf, ", get_mat_element_at(mat, l, c));
        printf("\n");
    }
    printf("]\n");
}

void free_mat_csr_t(mat_csr_t* mat) {
    free(mat->AA);
    free(mat->IA);
    free(mat->JA);
}

mat_csr_t mat_csr_t_from(double* mat_full, size_t row, size_t col) {
    mat_csr_t mat;
    mat.IA_l = row + 1;
    mat.IA = (size_t*)calloc(mat.IA_l, sizeof(size_t));
    mat.AA_l = 0;
    mat.JA = (size_t*)calloc(mat.AA_l, sizeof(size_t));
    mat.AA = (size_t*)calloc(mat.AA_l, sizeof(double));

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
        for(size_t i=r+1; i<mat.IA_l; i++) mat.IA[i] =  mat.AA_l;
    }
    return mat;
}

void print_vec(double* v, size_t n) {
    printf("[");
    for (size_t i = 0; i < n; i++) printf(i == n - 1 ? "%lf" : "%lf, ", v[i]);
    printf("]\n");
}

int main(int argc, char** argv) {
    double saad_example[25] = {1.0, 0, 0, 2, 0,  3,  4, 0, 5, 0, 6, 0, 7,
                               8,   9, 0, 0, 10, 11, 0, 0, 0, 0, 0, 12};
    double* saad_ptr = saad_example;
    mat_csr_t saad_csr = mat_csr_t_from(saad_ptr, 5, 5);
    print_mat_csr_t(&saad_csr);
    print_mat_csr_t_full(&saad_csr);

    free_mat_csr_t(&saad_csr);
}