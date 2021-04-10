#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort_vector.h>
#include <time.h>
#include <assert.h>
#include <gsl/gsl_linalg.h>

void matrix_print(gsl_matrix* A,FILE* fil) {
    for (int i = 0; i < (A->size1); i++) {
        for (int j = 0; j < (A->size2); j++) {
            double Aij = gsl_matrix_get(A, i, j);
            fprintf(fil, "%0.2f  ", Aij);
        }
        fprintf(fil, "\n");
    }
}

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
                        timesJ(A, p, q, theta);
                        Jtimes(A, p, q, -theta);
                        timesJ(V, p, q, theta);
                    }
            }
        }
    } while (changed!=0);
}


int main(){
    int n = 18;
    FILE* test = fopen("out.diagTest.txt","w");
    FILE* infWell = fopen("out.infWell.txt","w");
    FILE* infWellE = fopen("out.infWellEn.txt","w");

    gsl_matrix* A = gsl_matrix_alloc(n,n);
    gsl_matrix* V = gsl_matrix_alloc(n,n);
    gsl_matrix* Acpy = gsl_matrix_alloc(n,n);
    gsl_matrix* Vcpy = gsl_matrix_alloc(n,n);
    gsl_matrix* printMatrix = gsl_matrix_alloc(n,n);
    gsl_matrix* calc1 = gsl_matrix_alloc(n,n);
    gsl_matrix* calc2 = gsl_matrix_alloc(n,n);
    for(int i = 0; i<n; i++){
        for(int j=i; j<n; j++){
            double Aij = ((double) rand()/RAND_MAX*13);
            Aij -= ((double) rand()/RAND_MAX*13);
            gsl_matrix_set(A,i,j,Aij);
            gsl_matrix_set(A,j,i,Aij);
        }
    }
    gsl_matrix_transpose_memcpy(Acpy,A);
    gsl_matrix_set_identity(V);
    fprintf(test,"A=\n");
    matrix_print(A,test);


    jacobi_diag_opt(A,V);


    gsl_matrix_memcpy(Vcpy,V);
    fprintf(test,"D=\n");
    matrix_print(A,test);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.,Vcpy,V,0.,printMatrix);
    fprintf(test, "V^TV = 1?: \n");
    matrix_print(printMatrix,test);

    gsl_matrix_memcpy(printMatrix,A);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,Acpy,V,0.,calc1);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,V,calc1,-1.,printMatrix);
    fprintf(test, "V^T A V - D= 0 ?: \n");
    matrix_print(printMatrix,test);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.0,A,V,0.,calc1);
    gsl_matrix_memcpy(calc2,Acpy);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,V,calc1,-1.,calc2);
    fprintf(test, "V D V^T - A = 0?: \n");
    matrix_print(calc2,test);
    gsl_matrix_set_identity(calc1);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,calc1,A,-1.,Acpy);
    fprintf(test, " Anything happened? = D - A = 0?: \n");
    matrix_print(Acpy,test);
    gsl_matrix_free(A);
    gsl_matrix_free(Acpy);
    gsl_matrix_free(V);
    gsl_matrix_free(Vcpy);
    gsl_matrix_free(calc1);
    gsl_matrix_free(calc2);
    gsl_matrix_free(printMatrix);

    //HERE STARTS EXERCISE B
    int N = 20;
    double s = 1/(double) (N+1); //N+1 because we are never allowed to see the end, where there should be conditions
    gsl_matrix* H = gsl_matrix_calloc(N,N);
    gsl_matrix* W = gsl_matrix_alloc(N,N);
    gsl_matrix_set_identity(W);
    for(int i=0;i<N-1;i++){
        gsl_matrix_set(H,i,i,-2);
        gsl_matrix_set(H,i,i+1,1);
        gsl_matrix_set(H,i+1,i,1);
    }
    gsl_matrix_set(H,N-1,N-1,-2);
    gsl_matrix_scale(H,-1/s/s);
    jacobi_diag(H,W);

    gsl_vector_view v = gsl_matrix_diagonal(H);
    gsl_vector* w = gsl_vector_alloc(N);
    gsl_vector* indexvec = gsl_vector_alloc(N);
    gsl_vector_memcpy(w,&v.vector);
    for(int i =0; i<N; i++){
        gsl_vector_set(indexvec,i, i);
    }
    gsl_sort_vector2(w,indexvec);
    fprintf(infWellE,"Energies:\n");
    for (int k=0; k < n/3; k++){
        double exact = M_PI*M_PI*(k+1)*(k+1);
        double calculated = gsl_vector_get(w,k);
        fprintf(infWellE,"%i %g %g\n",k,calculated,exact);
    }

    fprintf(infWell,"%g %g %g %g %g %g %g\n", 0.,0.,0.,0.,0.,0.,0.);
    for(int i=0; i<N; i++){
        int j0 = gsl_vector_get(indexvec,0);
        int j1 = gsl_vector_get(indexvec,1);
        int j2 = gsl_vector_get(indexvec,2);
        fprintf(infWell,"%g %g %g %g %g %g %g\n", (i+1)*s,3.249*gsl_matrix_get(W,i,j0),3.249*gsl_matrix_get(W,i,j1),
                3.249*gsl_matrix_get(W,i,j2),sin(M_PI*(i+1)*s),sin(2*M_PI*(i+1)*s),-sin(3*M_PI*(i+1)*s));
    }
    fprintf(infWell,"%g %g %g %g %g %g %g\n", 1.,0.,0.,0.,0.,0.,0.);


    gsl_matrix_free(H);
    gsl_matrix_free(W);
    gsl_vector_free(w);
    gsl_vector_free(indexvec);

    //EXERCISE 3
    FILE* compare = fopen("out.compare.txt","w");
    double tiktok100 = 0;
    double gsltiktok100 = 0;
    double UTtiktok100 = 0;
    for(int i = -1; i<30; i++){
        int J =0;
        if(i==-1){
            J+=50;
        }
        else{
            J += 5*i+5;}
        gsl_matrix* M =gsl_matrix_alloc(J,J);
        gsl_matrix* G =gsl_matrix_alloc(J,J);
        gsl_matrix* Mut =gsl_matrix_alloc(J,J);
        gsl_matrix* Gut =gsl_matrix_alloc(J,J);
        gsl_matrix* gslM = gsl_matrix_alloc(J,J);
        gsl_matrix* gslV = gsl_matrix_alloc(J,J);
        gsl_vector* gslS = gsl_vector_alloc(J);
        for(int i = 0; i<J; i++){
        for(int j=i; j<J; j++){
            double Mij = ((double) rand()/RAND_MAX*13);
            Mij -= ((double) rand()/RAND_MAX*13);
            gsl_matrix_set(M,i,j,Mij);
            gsl_matrix_set(M,j,i,Mij);
        }
        }
        matrix_print(M,infWellE);
        gsl_matrix_memcpy(gslM,M);
        gsl_matrix_memcpy(Mut,M);
        gsl_matrix_set_identity(Gut);
        gsl_matrix_set_identity(G);
        gsl_matrix_set_identity(gslV);
        clock_t tik, tok, gsltik, gsltok, UTtik, UTtok;
        double timeUsed, gslTimeUsed, UTTimeUsed;
        tik=clock();
        jacobi_diag(M,G);
        tok=clock();

        UTtik=clock();
        jacobi_diag_UT(Mut,Gut);
        UTtok=clock();

        gsltik=clock();
        gsl_linalg_SV_decomp_jacobi(gslM, gslV, gslS);
        gsltok=clock();

        timeUsed = ((double) (tok - tik)) / CLOCKS_PER_SEC;
        gslTimeUsed = ((double) (gsltok - gsltik)) / CLOCKS_PER_SEC;
        UTTimeUsed = ((double) (UTtok - UTtik)) / CLOCKS_PER_SEC;
        if(i==-1){
            tiktok100+=timeUsed;
            gsltiktok100+=gslTimeUsed;
            UTtiktok100+=UTTimeUsed;
        }
        else {
            double N3 = pow(((double) J) / (50.), 3) * tiktok100;
            double gslN3 = pow(((double) J) / (50.), 3) * gsltiktok100;
            double UTN3 = pow(((double) J) / (50.), 3) * UTtiktok100;
            fprintf(compare, "%g %g %g %g %g %g %g\n", ((double) J), timeUsed, gslTimeUsed, UTTimeUsed, N3, gslN3,
                    UTN3);
        }
        gsl_matrix_free(M);
        gsl_matrix_free(G);
        gsl_matrix_free(Mut);
        gsl_matrix_free(Gut);
        gsl_matrix_free(gslM);
        gsl_matrix_free(gslV);
        gsl_vector_free(gslS);
        fclose(test);
        fclose(infWellE);
        fclose(infWell);
        fclose(compare);
    }
    return 0;
}