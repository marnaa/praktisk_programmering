#ifndef ROOTS_DECLARATIONS_H
#define ROOTS_DECLARATIONS_H
void matrix_print(gsl_matrix* A,FILE* fil);
void gs_decomp(gsl_matrix* A, gsl_matrix* R);
void backsub(gsl_matrix* R, gsl_vector* x);
void gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x);
void gs_invers(gsl_matrix* Q, gsl_matrix* R, gsl_matrix* X);
int bby_driver(void f(double x, gsl_vector* y, gsl_vector* dydx)
        ,double a, gsl_vector* ya, double b, double h,
               double acc, double eps);
void rkstep12(void f(double x, gsl_vector* y, gsl_vector* dydx), double x, gsl_vector* yx, double h, gsl_vector* yh, gsl_vector* err);

int adult_driver(void f(double x, gsl_vector* y, gsl_vector* dydx)
        ,double a, gsl_vector* ya, double **Yal,
                 int steps, double **Xal, double b, double h,
                 double acc, double eps);

#endif
