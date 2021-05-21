//
// Created by Martin on 21-05-2021.
//

#ifndef NEURAL_NETWORK_INTEGRATION_DECLARATIONS_H
#define NEURAL_NETWORK_INTEGRATION_DECLARATIONS_H
double integrater(double f(double x, ann* network), double a,
                  double b, ann* network,double acc, double eps, int* eta,double* error);
double CC_integrater(double f(double x, ann* network), double a,
                     double b, ann* network, double acc, double eps,int* eval, double* error);
#endif //NEURAL_NETWORK_INTEGRATION_DECLARATIONS_H
