//
// Created by Martin on 15-06-2021.
//
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef EXAM_ROUTINES_H

#define EXAM_ROUTINES_H

void lanczos(gsl_matrix* A, gsl_matrix* Q, gsl_matrix* H, gsl_vector* q0);

void matrix_print(gsl_matrix* A,FILE* fil);

void jacobi_diag_opt(gsl_matrix* A, gsl_matrix* V);

#endif //EXAM_ROUTINES_H
