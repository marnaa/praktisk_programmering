//
// Created by Martin on 11-05-2021.
//
#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#ifndef NEURAL_NETWORK_ANN_DECLARATIONS_H
#define NEURAL_NETWORK_ANN_DECLARATIONS_H
typedef struct {double(*f)(double); double(*f_diff)(double); double(*f_diffdiff)(double); double(*f_int)(double); gsl_vector* params; } ann;

double gaussWave(double x);

double gaussWave_diff(double x);
double gaussWave_diffdiff(double x);

double gaussWave_int(double x);

ann* ann_alloc(int n,double(*f)(double),double(*f_diff)(double),double(*f_diffdiff)(double),double(* f_int)(double ));

double ann_int(ann* network, double a, double b);

void   ann_free    (ann* network);

double ann_response(ann* network,double x);

double ann_diff(ann* network,double x);

double ann_diffdiff(ann* network,double x);

void ann_cost(ann* network, gsl_vector* xs, gsl_vector* ys);

void   ann_train   (ann* network,gsl_vector* xs,gsl_vector* ys);

double annWild_cost(ann* network, double diffeq_pow2(double responseofX),
double a, double b,double boundary_x,double boundary_y
,double boundary_ydot);

void   annWild_train  (ann* network, double diffeq_pow2(double responseofX, void* params),
                       double a, double b,double boundary_x,double boundary_y
        ,double boundary_ydot);
#endif //NEURAL_NETWORK_ANN_DECLARATIONS_H
