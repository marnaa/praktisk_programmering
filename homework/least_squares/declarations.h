#ifndef LEAST_SQUARES_DECLARATIONS_H
#define LEAST_SQUARES_DECLARATIONS_H
void matrix_print(gsl_matrix* A,FILE* fil);
void gs_decomp(gsl_matrix* A, gsl_matrix* R);
void backsub(gsl_matrix* R, gsl_vector* x);
void gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void gs_invers(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* X);

#endif //LEAST_SQUARES_DECLARATIONS_H
