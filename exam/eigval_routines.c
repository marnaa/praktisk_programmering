//
// Created by Martin on 16-06-2021.
//

#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <assert.h>

void timesJ(gsl_matrix* A, int p, int q, double theta){
    double c=cos(theta), s=sin(theta);
    for(int i=0; i<A->size1;i++){
        double AJ_ip = c*gsl_matrix_get(A,i,p)-s*gsl_matrix_get(A,i,q);
        double AJ_iq = c*gsl_matrix_get(A,i,q)+s*gsl_matrix_get(A,i,p);
        gsl_matrix_set(A,i,p,AJ_ip);
        gsl_matrix_set(A,i,q,AJ_iq);
    }
}
void Jtimes(gsl_matrix* A, int p, int q, double theta){
    double c=cos(theta), s=sin(theta);
    for(int i=0; i<A->size2;i++){
        double JA_pi = c*gsl_matrix_get(A,p,i)+s*gsl_matrix_get(A,q,i);
        double JA_qi = c*gsl_matrix_get(A,q,i)-s*gsl_matrix_get(A,p,i);
        gsl_matrix_set(A,p,i,JA_pi);
        gsl_matrix_set(A,q,i,JA_qi);
    }
}

//Jacobi diagonalisation, A is symmetric, Checking all elements
void jacobi_diag(gsl_matrix* A, gsl_matrix* V){
    int changed=0;
    int count=0;
    int n = A->size1;
    do{
        changed=0;
        count++;
        for(int p=0; p<n; p++) {
            for (int q = 0; q < n; q++) {
                if(p!=q){
                    double Apq = gsl_matrix_get(A, p, q);
                    double App = gsl_matrix_get(A, p, p);
                    double Aqq = gsl_matrix_get(A, q, q);
                    double theta = 0.5 * atan(2 * Apq/( Aqq - App));
                    double c = cos(theta), s = sin(theta);
                    double new_App = c * c * App - 2 * s * c * Apq + s * s * Aqq;
                    double new_Aqq = s * s * App + 2 * s * c * Apq + c * c * Aqq;
                    if (new_App!=App || new_Aqq!=Aqq) {
                        changed = 1;
                        timesJ(A, p, q, theta);
                        Jtimes(A, p, q, -theta);
                        timesJ(V, p, q, theta);
                    }
                }
                else{}
            }
        }
    } while (changed!=0);
}
void jacobi_diag_opt(gsl_matrix* A, gsl_matrix* V){
    int changed=0;
    int count=0;
    int n = A->size1;
    do{
        changed=0;
        count++;
        for(int p=0; p<n-1; p++) {
            for (int q = p+1; q < n; q++) {
                double Apq = gsl_matrix_get(A, p, q);
                double App = gsl_matrix_get(A, p, p);
                double Aqq = gsl_matrix_get(A, q, q);
                double theta = 0.5 * atan2(2 * Apq,( Aqq - App));
                double c = cos(theta), s = sin(theta);
                double new_App = c * c * App - 2 * s * c * Apq + s * s * Aqq;
                double new_Aqq = s * s * App + 2 * s * c * Apq + c * c * Aqq;
                if (new_App!=App || new_Aqq!=Aqq) {
                    changed = 1;
                    gsl_matrix_set(A,p,p,new_App);
                    gsl_matrix_set(A,q,q,new_Aqq);
                    gsl_matrix_set(A,p,q,0.0);
                    for(int i =0; i<p; i++){
                        double Aip = gsl_matrix_get(A,i,p);
                        double Aiq = gsl_matrix_get(A,i,q);
                        gsl_matrix_set(A,i,p,c*Aip-s*Aiq );
                        gsl_matrix_set(A,i,q,c*Aiq+s*Aip);
                    }
                    for(int i =p+1; i<q; i++){
                        double Api = gsl_matrix_get(A,p,i);
                        double Aiq = gsl_matrix_get(A,i,q);
                        gsl_matrix_set(A,p,i,c*Api-s*Aiq );
                        gsl_matrix_set(A,i,q,c*Aiq+s*Api);
                    }
                    for(int i =q+1; i<n; i++){
                        double Api = gsl_matrix_get(A,p,i);
                        double Aqi = gsl_matrix_get(A,q,i);
                        gsl_matrix_set(A,p,i,c*Api-s*Aqi);
                        gsl_matrix_set(A,q,i,c*Aqi+s*Api);
                    }
                    for(int i=0; i<n; i++){
                        double Vip=gsl_matrix_get(V,i,p);
                        double Viq = gsl_matrix_get(V,i,q);
                        gsl_matrix_set(V,i,p,c*Vip-s*Viq);
                        gsl_matrix_set(V,i,q,c*Viq+s*Vip);
                    }
                }
            }
        }
    } while (changed!=0);
}
