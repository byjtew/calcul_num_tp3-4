#include "lib_poisson1D.h"

int nbpoints, la, info;
int ku, kl, kv, lab;
int *row_ipiv, *col_ipiv;
int row_NRHS, col_NRHS;
double T0, T1;
double *row_RHS, *col_RHS, *ROW_EX_SOL, *COL_EX_SOL, *X;
double *AB;
double temp, relres;

void initialize() {
    row_NRHS = 1;
    col_NRHS = 1;
    nbpoints = 12;
    la = nbpoints - 2;
    T0 = -5.0;
    T1 = 5.0;

    row_RHS = (double *)malloc(sizeof(double) * la);
    col_RHS = (double *)malloc(sizeof(double) * la);
    ROW_EX_SOL = (double *)malloc(sizeof(double) * la);
    COL_EX_SOL = (double *)malloc(sizeof(double) * la);
    X = (double *)malloc(sizeof(double) * la);

    set_grid_points_1D(X, &la);
    set_dense_RHS_DBC_1D(row_RHS, &la, &T0, &T1);
    set_dense_RHS_DBC_1D(col_RHS, &la, &T0, &T1);
    set_analytical_solution_DBC_1D(ROW_EX_SOL, X, &la, &T0, &T1);
    set_analytical_solution_DBC_1D(COL_EX_SOL, X, &la, &T0, &T1);

    kv = 1;
    ku = 1;
    kl = 1;
    lab = kv + kl + ku + 1;

    AB = (double *)malloc(sizeof(double) * lab * la);

    /* working array for pivot used by LU Factorization */
    row_ipiv = (int *)calloc(la, sizeof(int));
    col_ipiv = (int *)calloc(la, sizeof(int));
}

int main(int argc, char *argv[]) {
    printf("--- Poisson 1D: Difference test between ROW & COL MAJOR ---\n\n");

    initialize();

    /* 
        ROW SECTION 
    */

    set_GB_operator_rowMajor_poisson1D(AB, &lab, &la, &kv);
    // write_GB_operator_rowMajor_poisson1D(AB, &lab, &la, "AB_row.dat");

    info = LAPACKE_dgbsv(LAPACK_ROW_MAJOR, la, kl, ku, row_NRHS, AB, la,
                         col_ipiv, row_RHS, row_NRHS);

    /* Relative residual */
    temp = cblas_ddot(la, row_RHS, 1, row_RHS, 1);
    temp = sqrt(temp);
    cblas_daxpy(la, -1.0, row_RHS, 1, ROW_EX_SOL, 1);
    relres = cblas_ddot(la, ROW_EX_SOL, 1, ROW_EX_SOL, 1);
    relres = sqrt(relres);
    relres = relres / temp;

    printf("\nThe relative residual error for ROW_MAJOR is relres = %e\n",
           relres);

    /* =================== */

    // Reset AB
    for (int i = 0; i < la * lab; i++) AB[i] = .0;

    /* 
        COL SECTION 
    */

    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
    // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_col.dat");

    info = LAPACKE_dgbsv(LAPACK_COL_MAJOR, la, kl, ku, col_NRHS, AB, lab,
                         row_ipiv, col_RHS, la);

    printf("\n INFO DGBSV = %d\n", info);

    double cum_diff = .0;
    for (int i = 0; i < la; i++) cum_diff += fabs(row_RHS[i] - col_RHS[i]);

    printf("Difference between row result and col result: %lf\n", cum_diff);

    /* Relative residual */
    temp = cblas_ddot(la, row_RHS, 1, row_RHS, 1);
    temp = sqrt(temp);
    cblas_daxpy(la, -1.0, row_RHS, 1, COL_EX_SOL, 1);
    relres = cblas_ddot(la, COL_EX_SOL, 1, COL_EX_SOL, 1);
    relres = sqrt(relres);
    relres = relres / temp;

    printf("\nThe relative residual error for COL_MAJOR is relres = %e\n",
           relres);

    /* =================== */



    free(row_RHS);
    free(col_RHS);
    free(ROW_EX_SOL);
    free(COL_EX_SOL);
    free(X);
    free(AB);
    free(row_ipiv);
    free(col_ipiv);

    printf("\n\n--- End of the test ---\n");
}