#include "routines.h"
#include <gsl/gsl_sort_vector.h>

#define RAND (double) rand()/RAND_MAX
int main(){

    //Defining file for eigenvalue check
    FILE* eigval_check = fopen("out.eigval_check.txt","w");

    int n = 10;
    int m = 6;

    int number_size = 10;

    // Defining the matrices needed for the computation
    gsl_matrix* A = gsl_matrix_calloc(n, n);
    gsl_matrix* H = gsl_matrix_calloc(n, n);
    gsl_matrix* Q = gsl_matrix_calloc(n, n);
    gsl_matrix* QH = gsl_matrix_calloc(n, n);
    gsl_matrix* cpyA = gsl_matrix_calloc(n,n);

    //Defining matrices for a Kyrlow space with m<n
    gsl_matrix* Hm = gsl_matrix_calloc(m, m);
    gsl_matrix* Qm = gsl_matrix_calloc(n, m);
    gsl_matrix* QHm = gsl_matrix_calloc(n, m);
    gsl_matrix* a = gsl_matrix_calloc(n,n);



    //Defining initial guess
    gsl_vector* q0 = gsl_vector_calloc(n);

    //Making A a symmetric matrix
    for(int i=0; i<n; i++){
        for(int j=i; j<n; j++){
            double Aij = number_size*RAND;
            Aij -= number_size*RAND;
            gsl_matrix_set(A,i,j,Aij);
            gsl_matrix_set(A,j,i,Aij);
        }
    }
    gsl_matrix_memcpy(cpyA,A);
    gsl_matrix_memcpy(a,A);
    gsl_vector_set(q0, 1,1.);

    //Calling lanczos
    lanczos(A,Q,H,q0);
    lanczos(A,Qm,Hm,q0);
    printf("TESTING THAT THE MATRICES ARE FULFILLING THE RELATIONS GIVEN BY THEORY\n\n");

    printf("m=n=10\n");
    printf("A:\n");
    matrix_print(A,stdout);
    printf("\n");
    printf("T:\n");
    matrix_print(H,stdout);
    printf("\n");
    printf("V:\n");
    matrix_print(Q,stdout);
    printf("\n");
    // Making the proper tests

    printf("V T V^T - A = 0 ? (only valid m=n): \n");
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,Q,H,0.,QH);
    gsl_blas_dgemm(CblasNoTrans,CblasTrans,1.,QH,Q,-1.,A);
    matrix_print(A,stdout);
    printf("\n");


    printf("V^T A V - T = 0 ?:\n");
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,cpyA,Q,0.,QH);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,Q,QH,-1.,H);
    matrix_print(H,stdout);
    printf("\n");

    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,Q, Q,0.,H);
    printf("V^TV=1 ?:\n");
    matrix_print(H,stdout);
    printf("\n\n");

    printf("m=6, n=10\n");

    printf("A:\n");
    matrix_print(cpyA,stdout);
    printf("\n");
    printf("T:\n");
    matrix_print(Hm,stdout);
    printf("\n");
    printf("V:\n");
    matrix_print(Qm,stdout);
    printf("\n");

    printf("V^T A V - T = 0 ?:\n");
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,cpyA,Qm,0.,QHm);
    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,Qm,QHm,-1.,Hm);
    matrix_print(Hm,stdout);
    printf("\n");

    gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.,Qm, Qm,0.,Hm);
    printf("V^TV=1 ?:\n");
    matrix_print(Hm,stdout);
    printf("\n\n");


    fprintf(eigval_check,"CHECKING HOW WELL THE EIGENVALUES OF T APPROXIMATES THE HIGHEST EIGENVALUE(S) OF A\n"
                         "WITH DIFFERENT SIZES OF KYRLOV BASES [RANK(A)=10] \n\n");

    //Finding the A eigenvalues
    gsl_matrix* va = gsl_matrix_calloc(n,n);
    gsl_matrix* acpy = gsl_matrix_calloc(n,n);
    gsl_matrix_memcpy(acpy,a);
    jacobi_diag_opt(acpy,va);

    //Sorting the eigenvalues of A after size
    gsl_vector_view v = gsl_matrix_diagonal(acpy);
    gsl_vector* p = gsl_vector_alloc(n);
    gsl_vector* p_index = gsl_vector_alloc(n);

    gsl_vector_memcpy(p,&v.vector);
    for(int j =0; j<n; j++){
        gsl_vector_set(p_index,j, j);
    }
    gsl_sort_vector2(p,p_index);

    for( int i = 1; i<=n; i++){
        //Defining matrices for eigenvalue comparison
        gsl_matrix* h = gsl_matrix_calloc(i,i);
        gsl_matrix* q = gsl_matrix_calloc(n,i);
        gsl_matrix* vh = gsl_matrix_calloc(i,i);
        gsl_vector* w = gsl_vector_alloc(i);
        gsl_vector* indexvec = gsl_vector_alloc(i);
        //Doing lanczos and finding eigenvalues of H through the jacobi routines from homework
        lanczos(a,q,h,q0);
        jacobi_diag_opt(h,vh);

        //Sorting the eigenvalues after size
        fprintf(eigval_check,"The eigenvalue at a subspace of %i:\n",i);
        gsl_vector_view v = gsl_matrix_diagonal(h);


        gsl_vector_memcpy(w,&v.vector);
        for(int j =0; j<i; j++){
            gsl_vector_set(indexvec,j, j);
        }
        gsl_sort_vector2(w,indexvec);

        //
        //Printing and comparing the eigenvalues of h and A
        fprintf(eigval_check,"Eigenvalues (T, A):\n");
        for(int k=0; k<i; k++){
            double h_eig = gsl_vector_get(w,i-(k+1));
            double a_eig = gsl_vector_get(p,n-(k+1));
            fprintf(eigval_check,"%i %g %g\n",(k+1),h_eig,a_eig);

        }


        gsl_matrix_free(h);
        gsl_matrix_free(q);
        gsl_matrix_free(vh);
        gsl_vector_free(w);
        fprintf(eigval_check," \n");
        fprintf(eigval_check," \n");
    }
    //Freeing and closing hopefully everything
    gsl_matrix_free(va);
    gsl_matrix_free(acpy);
    gsl_matrix_free(A);
    gsl_matrix_free(H);
    gsl_matrix_free(cpyA);
    gsl_matrix_free(QH);
    gsl_matrix_free(Q);
    gsl_matrix_free(Hm);
    gsl_matrix_free(a);
    gsl_matrix_free(QHm);
    gsl_matrix_free(Qm);
    gsl_vector_free(p);
    gsl_vector_free(p_index);
    fclose(eigval_check);

    return 0;
}
