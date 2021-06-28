#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <assert.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>



void matrix_print(gsl_matrix* A,FILE* fil){
    for(int i=0; i<(A->size1);i++){
        for(int j=0; j<(A->size2); j++){
            double Aij = gsl_matrix_get(A,i,j);
            fprintf(fil, "%0.3f  ",Aij);
        }
        fprintf(fil,"\n");
    }
}
//Gram-Schmidt orthogonalisation function
void gs_decomp(gsl_matrix* A, gsl_matrix* R) {
    gsl_vector *ai = gsl_vector_alloc((A->size1));
    gsl_vector *qi = gsl_vector_alloc((A->size1));
    //Putting normalised ai's into Ri
    for (int i = 0; i < (A->size2); i++) {
        //finding norm of a_i
        gsl_matrix_get_col(ai, A, i);
        gsl_vector_memcpy(qi,ai);
        double ai_norm = gsl_blas_dnrm2(ai);
        assert(ai_norm>0);
        //Setting R_ii = norm(a_i, a_i)
        gsl_matrix_set(R, i, i, ai_norm);
        gsl_vector_scale(qi,1/ai_norm);
        gsl_matrix_set_col(A,i,qi);
        for(int j=i+1; j<(A->size2); j++){
            //a_i==a_j
            gsl_matrix_get_col(ai,A,j);
            double Rij=0;
            //Rij=q^T a_i
            gsl_blas_ddot(qi,ai,&Rij);
            //a_j=-q_i+a_i
            gsl_blas_daxpy(-Rij,qi,ai);
            //set a_j coloumn
            gsl_matrix_set(R,i,j,Rij);
            gsl_matrix_set(R,j,i,0.);
            gsl_matrix_set_col(A,j,ai);
        }

    }
    gsl_vector_free(ai);
    gsl_vector_free(qi);
}

void backsub(gsl_matrix* R, gsl_vector* x){
    for(int i=(x->size)-1; i>=0;i--){
        double QTbi = gsl_vector_get(x,i);
        double Rii = gsl_matrix_get(R,i,i);
        double xi = QTbi;
        for(int j=i+1; j<x->size;j++){
            double Rij = gsl_matrix_get(R,i,j);
            double xj = gsl_vector_get(x,j);
            xi-=Rij*xj;
        }
        xi /= Rii;
        gsl_vector_set(x,i,xi);
    }
}

void gs_solve(gsl_matrix* Q, gsl_matrix* R, gsl_vector* b, gsl_vector* x){
    gsl_blas_dgemv(CblasTrans,1.0,Q,b,0.0,x);
    backsub(R,x);
}

void gs_invers(gsl_matrix *Q, gsl_matrix* R, gsl_matrix* X){
    //Need to solve QRB=1, in principle also BQR=1, but not nice form
    gsl_matrix_transpose_memcpy(X,Q);
    for(int i=0; i<(Q->size2); i++){
        gsl_vector_view qi = gsl_matrix_column(X,i);
        backsub(R, &qi.vector);
    }
}


int main(){
    FILE* gstest = fopen("out.gstest.txt","w");
    int n = 10;
    int tall = 13;
    int m = 10;
    gsl_matrix *A = gsl_matrix_alloc(n,m);
    gsl_matrix *Atall = gsl_matrix_alloc(tall,m);
    gsl_matrix *Rtall = gsl_matrix_alloc(m,m);
    gsl_matrix *Atall_test = gsl_matrix_calloc(m,m);
    gsl_matrix *Atall_copy = gsl_matrix_alloc(tall,m);


    gsl_matrix *R = gsl_matrix_alloc(m,m);
    gsl_matrix *X = gsl_matrix_alloc(m,m);
    gsl_vector *x =gsl_vector_alloc(m);
    gsl_vector *b =gsl_vector_alloc(m);

    //Setting values for the tall matrix
    for(int j=0; j<(Atall->size2); j++){
        for(int i = 0; i<(Atall->size1);i++){
            double Aij = (double) rand()/RAND_MAX*40;
            Aij-=(double) rand()/RAND_MAX*40;
            gsl_matrix_set(Atall,i,j,Aij);
        }
    }

    //Setting values for the square matrix
    for(int j=0; j<(A->size2); j++){
        double bj = (double) rand()/RAND_MAX*40;
        bj-=(double) rand()/RAND_MAX*40;
        gsl_vector_set(b,j,bj);
        for(int i = 0; i<(A->size1);i++){
            double Aij = (double) rand()/RAND_MAX*40;
            Aij-=(double) rand()/RAND_MAX*40;
            gsl_matrix_set(A,i,j,Aij);
        }
    }

    //Copy matrices for checking the proper relations
    gsl_vector *bcopy = gsl_vector_alloc(m);
    gsl_matrix *AA = gsl_matrix_alloc(n,m);
    gsl_matrix *XX = gsl_matrix_alloc(n,m);
    gsl_matrix *XXX = gsl_matrix_alloc(n,m);
    gsl_matrix *Acopy = gsl_matrix_alloc(n,m);
    gsl_matrix *Qtest = gsl_matrix_alloc(m,m);
    gsl_matrix_memcpy(Acopy,A);
    gsl_vector_memcpy(bcopy,b);
    gsl_matrix_memcpy(AA,A);
    gsl_matrix_memcpy(Atall_copy,Atall);

    //The actual decomposition
    //tall matrix
    gs_decomp(Atall,Rtall);

    //square matrix
    gs_decomp(A,R);
    gs_solve(A,R,b,x);
    gs_invers(A,R,X);


    //Making the tall test Q^T Q =1 QR=A, R upper triangle
    fprintf(gstest,"TALL TEST\n \n ");
    fprintf(gstest,"R upper triangle ?\n ");
    matrix_print(R,gstest);

    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,A,A,0.,Atall_test);

    fprintf(gstest,"Q^T Q = 1 ? \n");
    matrix_print(Atall_test,gstest);

    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,Atall,Rtall,-1.,Atall_copy);

    fprintf(gstest,"QR-A = 0 ? \n");
    matrix_print(Atall_copy,gstest);
    fprintf(gstest,"\n");


    //Making the square tests

    fprintf(gstest,"SQUARE/INVERSE/SOLUTION TEST\n \n ");
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,A,R,-1.,Acopy);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,A,A,0.,Qtest);
    gsl_blas_dgemv(CblasNoTrans,1.0,AA,x,-1.,bcopy);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,AA,X,0.,XX);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,X,AA,0.,XXX);

    fprintf(gstest,"A:\n");
    matrix_print(A,gstest);

    fprintf(gstest,"QR-A=0?:\n");
    matrix_print(Acopy,gstest);

    fprintf(gstest,"(Q^T)Q=1?:\n");
    matrix_print(Qtest,gstest);

    fprintf(gstest,"R uppertriangle ?:\n");
    matrix_print(R,gstest);

    fprintf(gstest,"b:\n");
    gsl_vector_fprintf(gstest,b,"%0.3g");

    fprintf(gstest,"Ax-b=0 ?:\n");
    gsl_vector_fprintf(gstest,bcopy,"%0.3g");

    fprintf(gstest,"A(A^(-1))=1 ?:\n");
    matrix_print(XX,gstest);

    fprintf(gstest,"(A^(-1))A=1 ?:\n");
    matrix_print(XXX,gstest);

    FILE* time = fopen("out.time.txt","w");
    FILE* compare = fopen("out.compare.txt","w");
    int count = 0;
    double tiktok100 = 0;
    while(count<100){
        int N=0;
        if(count==0){
            N += 100;
        }
        else {
            N += (double) rand() / RAND_MAX * 500;
        }
        assert(N<501);
        gsl_matrix* M =gsl_matrix_alloc(N,N);
        gsl_matrix* G =gsl_matrix_alloc(N,N);
        gsl_matrix* gslM = gsl_matrix_alloc(N,N);
        gsl_vector* gslTau = gsl_vector_alloc(N);
        for(int j=0; j<(M->size2); j++){
            for(int i = 0; i<(M->size1);i++){
                double Mij = (double) rand()/RAND_MAX*40;
                Mij-=(double) rand()/RAND_MAX*40;
                gsl_matrix_set(M,i,j,Mij);
            }
        }
        gsl_matrix_memcpy(gslM,M);
        clock_t tik, tok, gsltik, gsltok;
        double timeUsed, gslTimeUsed;
        tik=clock();
            gs_decomp(M,G);
        tok=clock();

        gsltik=clock();
            gsl_linalg_QR_decomp(gslM,gslTau);
        gsltok=clock();

        timeUsed = ((double) (tok - tik)) / CLOCKS_PER_SEC;
        gslTimeUsed = ((double) (gsltok - gsltik)) / CLOCKS_PER_SEC;
        if(count==0){
            tiktok100+=timeUsed;
        }
        double N3 = pow(((double) N)/(100),3)*tiktok100;
        fprintf(time,"%g %g %g\n",((double)N),timeUsed,N3);
        fprintf(compare,"%g %g %g\n",((double)N),timeUsed,gslTimeUsed);
        count++;
        gsl_matrix_free(M);
        gsl_matrix_free(G);
        gsl_matrix_free(gslM);
        gsl_vector_free(gslTau);
    }
    //Freeing hopefully everything
    gsl_vector_free(bcopy);
    gsl_matrix_free(AA);
    gsl_matrix_free(XX);
    gsl_matrix_free(XXX);
    gsl_matrix_free(Acopy);
    gsl_matrix_free(Qtest);
    gsl_matrix_free(A);
    gsl_matrix_free(Atall);
    gsl_matrix_free(Rtall);
    gsl_matrix_free(Atall_test);
    gsl_matrix_free(Atall_copy);
    gsl_matrix_free(R);
    gsl_matrix_free(X);
    gsl_vector_free(x);
    gsl_vector_free(b);

    fclose(gstest);




    return 0;
}