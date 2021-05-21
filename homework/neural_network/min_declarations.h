#include <math.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#ifndef NEURAL_NETWORK_DECLARATIONS_H
#define NEURAL_NETWORK_DECLARATIONS_H

typedef struct {double(*f)(double); double(*f_diff)(double); double(*f_diffdiff)(double); double(*f_int)(double); gsl_vector* params; } ann;

void ann_amoeba( double cost(ann* network, gsl_vector* xs, gsl_vector* ys),
                 ann* network, gsl_vector* xs, gsl_vector* ys, double eps);
void annWild_amoeba( double cost(ann* network, double diffeq_pow2(double responseofX, ann* network),
                                 double a, double b,double boundary_x,double boundary_y
        ,double boundary_ydot), ann* network, double diffeq_pow2(double responseofX, ann* network),
                     double a, double b,double boundary_x,double boundary_y
        ,double boundary_ydot, double eps);
void amoeba( double f(gsl_vector* x),
             gsl_vector* x, gsl_vector* step, double eps);
void qnewton( double f(gsl_vector* x),
              gsl_vector* x, double eps);

#endif //NEURAL_NETWORK_DECLARATIONS_H
