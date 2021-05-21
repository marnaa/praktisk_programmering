#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<float.h>
#include<math.h>
static int N_MAX;
void matrix_print(gsl_matrix* A,FILE* fil){
    for(int i=0; i<(A->size1);i++){
        for(int j=0; j<(A->size2); j++){
            double Aij = gsl_matrix_get(A,i,j);
            fprintf(fil, "%0.3f  ",Aij);
        }
        fprintf(fil,"\n");
    }
}

void gradient(double f(gsl_vector* x),gsl_vector* x, gsl_vector* gradient,
              gsl_vector* Dx){
    int n = x -> size;
    gsl_vector_memcpy(Dx,x);
    for (int j = 0; j < n; j++) {
        double xj = gsl_vector_get(x, j);
        gsl_vector_set(Dx, j, xj+sqrt(DBL_EPSILON));
        double J_j = (f(Dx)-f(x))/sqrt(DBL_EPSILON);
        gsl_vector_set(gradient, j, J_j);
        gsl_vector_memcpy(Dx,x);
    }
}
void qnewton( double f(gsl_vector* x),
              gsl_vector* x, double eps){
    int n = x->size;
    double alpha = 1e-4;
    double prikprod;
    double prikprodsy;
    gsl_matrix* B = gsl_matrix_alloc(n,n);
    gsl_vector* grad = gsl_vector_alloc(n);
    gsl_vector* xs = gsl_vector_alloc(n);
    gsl_vector* grad_vec = gsl_vector_alloc(n);
    gsl_vector* s = gsl_vector_alloc(n);
    gsl_vector* grad_xs=gsl_vector_alloc(n);
    gsl_vector* y=gsl_vector_alloc(n);
    gsl_vector* u=gsl_vector_alloc(n);
    N_MAX = 0;
    gsl_matrix_set_identity(B);
    gradient(f,x,grad,grad_vec);
    while(gsl_blas_dnrm2(grad)>eps && N_MAX<10000){
        N_MAX++;
        double lambda = 1;
        gsl_blas_dgemv(CblasNoTrans,-1.,B,grad,0.,s);
        for(int i=0; i<n; i++){
            gsl_vector_set(xs,i, gsl_vector_get(s,i)+gsl_vector_get(x,i));
        }
        gsl_vector_memcpy(xs,x);
        gsl_vector_add(xs,s);
        gsl_blas_ddot(s,grad,&prikprod);
        double fx = f(x);
        while(f(xs)>fx+alpha*prikprod){
            if(lambda<sqrt(DBL_EPSILON)){
                gsl_matrix_set_identity(B);
                break;
            }
            lambda/=2;
            //printf("%g %g\n",f(xs), f(x)+alpha*prikprod);
            gsl_vector_scale(s,0.5);
            gsl_vector_memcpy(xs,x);
            gsl_vector_add(xs,s);
            gsl_blas_ddot(s,grad,&prikprod);
        }
        gradient(f,xs,grad_xs,grad_vec);
        gsl_vector_memcpy(y,grad_xs);
        gsl_blas_daxpy(-1,grad,y); //y = grad(x+s)-grad(x)
        gsl_vector_memcpy(u,s);
        gsl_blas_dgemv(CblasNoTrans,-1.,B,y,1.,u);
        gsl_blas_ddot(s,y,&prikprodsy);
        if(fabs(prikprodsy)>1e-12){
            gsl_blas_dger(1./prikprodsy,u,s,B);
        }
        gsl_vector_memcpy(x,xs);
        gsl_vector_memcpy(grad ,grad_xs);
    }
    gsl_vector_free(grad);
    gsl_matrix_free(B);
    gsl_vector_free(grad_vec);
    gsl_vector_free(xs);
    gsl_vector_free(y);
    gsl_vector_free(grad_xs);
    gsl_vector_free(u);
    gsl_vector_free(s);
}

/*
 * Implementing amoeba minimization
 * so
 * that's it
 */
void hiLoCent(gsl_matrix* simplex, gsl_vector* F_val, gsl_vector* centroid, int* hi, int* lo){
    //finding highest and lowest value
    int n = centroid -> size;
    *hi = 0;
    *lo = 0;
    double highest = gsl_vector_get(F_val,*hi);
    double lowest = gsl_vector_get(F_val,*lo);
    for(int i = 1; i<n+1;i++){
        double fi = gsl_vector_get(F_val,i);
        if(fi>highest){
            *hi = i;
            highest = fi;
        }
        else if(fi<lowest){
            *lo = i;
            lowest = fi;
        }
    }
    //finding centroid
    for(int i =0; i<n; i++){
        double sum = 0;
        for(int j =0; j<n+1;j++){
            if(j != *hi){
                sum+=gsl_matrix_get(simplex,i,j)/n;
            }
            gsl_vector_set(centroid,i,sum);
        }
    }
}

//Defining ann here so to make the adjustments needed for ann.c later on
typedef struct { double(*f)(double); double(*f_diff)(double); double(*f_diffdiff)(double); double(*f_int)(double); gsl_vector* params; } ann;

void ann_init_vals(double cost(ann* network, gsl_vector* xs, gsl_vector* ys), gsl_matrix* simplex, ann* init,
               gsl_vector* xs, gsl_vector* ys, gsl_vector* F_val, gsl_vector* centroid, int* hi, int* lo){
    int n = centroid -> size;
    for(int i = 0; i<n+1; i++){
        gsl_vector_view x = gsl_matrix_column(simplex,i);
        init->params = &x.vector;
        /*
        printf("x-vector in params?:\n");
        gsl_vector_fprintf(stdout,(init->params),"%g");
        printf("cost function:\n");
        printf("%g\n",cost(init,xs,ys));
        */
        gsl_vector_set(F_val,i,cost(init,xs,ys));
    }
    hiLoCent(simplex,F_val,centroid,hi,lo);
}

double size(gsl_matrix* simplex, int lo){
    //finding the greatest distance between the points in the simplex
    double s=0.;
    double distance;
    int n = simplex -> size1;
    gsl_vector* dist = gsl_vector_alloc(n);
    gsl_vector_view lowest = gsl_matrix_column(simplex,lo);
    for(int i =0; i<n+1;i++){
        gsl_vector_view xi = gsl_matrix_column(simplex,i);
        gsl_vector_memcpy(dist,&lowest.vector);
        gsl_vector_sub(dist,&xi.vector);
        distance = gsl_blas_dnrm2(dist);
        if(distance>s){
            s=distance;
        }
    }
    gsl_vector_free(dist);
    return s;
}

void reflect(gsl_vector* highest, gsl_vector* centroid, gsl_vector* reflected){
    gsl_vector_memcpy(reflected,highest);
    gsl_vector_scale(reflected,-1.);
    gsl_blas_daxpy(2.,centroid,reflected);
}
void expand(gsl_vector* highest, gsl_vector* centroid, gsl_vector* expanded){
    gsl_vector_memcpy(expanded,highest);
    gsl_vector_scale(expanded,-2.);
    gsl_blas_daxpy(3.,centroid,expanded);
}
void contract(gsl_vector* highest, gsl_vector* centroid, gsl_vector* contracted){
    gsl_vector_memcpy(contracted,highest);
    gsl_vector_scale(contracted,0.5);
    gsl_blas_daxpy(0.5,centroid,contracted);
}

void reduce(gsl_matrix* simplex, int lo){
    int n = simplex -> size1;
    for(int i=0; i<n+1;i++){
        if(i!=lo){
            gsl_vector_view lowest = gsl_matrix_column(simplex,lo);
            gsl_vector_view i_vec = gsl_matrix_column(simplex,i);
            gsl_vector_add(&i_vec.vector,&lowest.vector);
            gsl_vector_scale(&i_vec.vector,0.5);
        }
    }
}

void ann_amoeba( double cost(ann* network, gsl_vector* xs, gsl_vector* ys),
             ann* network, gsl_vector* xs, gsl_vector* ys, double eps){
    int n = (network->params)->size;
    gsl_vector* step = gsl_vector_calloc(n);
    gsl_vector_memcpy(step,(network->params));
    for(int i=0; i<n; i++){
        double stepi = gsl_vector_get(step,i);
        gsl_vector_set(step,i,0.7*stepi);
    }
    N_MAX = 0;
    gsl_matrix* simplex = gsl_matrix_alloc(n,n+1);
    gsl_vector* centroid = gsl_vector_alloc(n);

    gsl_vector* init_vec = gsl_vector_alloc(n);
    ann* init = malloc(sizeof(ann));
    init->f = network->f;
    init->params = init_vec;

    ann* ann_p1 = malloc(sizeof(ann));
    gsl_vector* p1 = gsl_vector_alloc(n);
    ann_p1->f = (network->f);
    ann_p1->params = p1;

    ann* ann_p2 = malloc(sizeof(ann));
    gsl_vector* p2 = gsl_vector_alloc(n);
    ann_p2->f = (network->f);
    ann_p2 -> params = p2;

    gsl_vector* F_val = gsl_vector_alloc(n+1);
    int hi, lo;

    //starting point
    for(int i = 0; i<n+1; i++){
        gsl_matrix_set_col(simplex,i,(network->params));
    }
    for(int i = 0; i<n; i++){
        double ii = gsl_matrix_get(simplex,i,i);
        double stepi = gsl_vector_get(step,i);
        gsl_matrix_set(simplex,i,i, ii+stepi);
    }

    //finding high, low, centroid
    ann_init_vals(cost,simplex,init, xs, ys,F_val,centroid,&hi,&lo);
    //printf("%g\n",size(simplex,lo));
    //matrix_print(simplex,stdout);
    //printf("%g\n",size(simplex,lo));
    while(size(simplex,lo)>eps && N_MAX<1e7){
        N_MAX++;
        hiLoCent(simplex,F_val,centroid,&hi,&lo);
        gsl_vector_view highest = gsl_matrix_column(simplex,hi);
        /*
        printf("N_MAX: %i \n ",N_MAX);

        printf("centroid: \n ");
        gsl_vector_fprintf(stdout,centroid,"%g");
        printf("Fvals: \n ");
        gsl_vector_fprintf(stdout,F_val,"%g");

        printf("Simplex: \n ");
        matrix_print(simplex,stdout);
        printf("highest = %i , lowest = %i\n", hi, lo);
        printf("Highest: \n");
        gsl_vector_fprintf(stdout,&highest.vector,"%g");
        */

        reflect(&highest.vector, centroid, p1);


        double f_re = cost(ann_p1,xs,ys);
        //printf("P1: \n");
        //gsl_vector_fprintf(stdout, p1 ,"%g");

        if(f_re< gsl_vector_get(F_val,lo)){
            expand(&highest.vector,centroid,p2);
            double f_ex = cost(ann_p2,xs,ys);
            if(f_ex<f_re){
                gsl_vector_memcpy(&highest.vector,p2);
                gsl_vector_set(F_val,hi,f_ex);
            }
            else{
                gsl_vector_memcpy(&highest.vector,p1);
                gsl_vector_set(F_val,hi,f_re);
            }
        }
        else{
            if(f_re<gsl_vector_get(F_val,hi)){
                gsl_vector_memcpy(&highest.vector,p1);
                gsl_vector_set(F_val,hi,f_re);
            }

            else {
                contract(&highest.vector, centroid, p1);
                double f_co = cost(ann_p1,xs,ys);
                if (f_co < gsl_vector_get(F_val, hi)) {
                    gsl_vector_memcpy(&highest.vector, p1);
                    gsl_vector_set(F_val, hi, f_co);
                }
                else {
                    reduce(simplex, lo);
                    ann_init_vals(cost,simplex,init, xs, ys,F_val,centroid,&hi,&lo);
                }
            }

        }
    }
    gsl_vector_view final_low = gsl_matrix_column(simplex,lo);
    //printf("final\n");
    //gsl_vector_fprintf(stdout,&final_low.vector,"%g");
    gsl_vector_memcpy((network->params),&final_low.vector);
    //printf("params\n");
    //gsl_vector_fprintf(stdout,(network->params),"%g");
    printf("N_MAX = %i\n",N_MAX);
    gsl_vector_free(centroid);
    gsl_vector_free(F_val);
    gsl_vector_free(p1);
    gsl_vector_free(p2);
    free(ann_p1);
    free(ann_p2);
    free(init);
    gsl_matrix_free(simplex);
}

//FOR UNSUPERVISED TRAINING
void annWild_init_vals(double cost(ann* network, double diffeq_pow2(double responseofX, ann* network),
                                   double a, double b,double boundary_x,double boundary_y
        ,double boundary_ydot), gsl_matrix* simplex, ann* init, double diffeq_pow2(double responseofX, ann* network),
                       double a, double b,double boundary_x,double boundary_y
        ,double boundary_ydot, gsl_vector* F_val, gsl_vector* centroid, int* hi, int* lo){
    int n = centroid -> size;
    for(int i = 0; i<n+1; i++){
        gsl_vector_view x = gsl_matrix_column(simplex,i);
        init->params = &x.vector;
        /*
        printf("x-vector in params?:\n");
        gsl_vector_fprintf(stdout,(init->params),"%g");
        printf("cost function:\n");
        printf("%g\n",cost(init,xs,ys));
        */
        gsl_vector_set(F_val,i,cost(init, diffeq_pow2,
        a, b, boundary_x, boundary_y
        ,boundary_ydot));
    }
    hiLoCent(simplex,F_val,centroid,hi,lo);
}
void annWild_amoeba( double cost(ann* network, double diffeq_pow2(double responseofX, ann* network),
                                        double a, double b,double boundary_x,double boundary_y
        ,double boundary_ydot), ann* network, double diffeq_pow2(double responseofX, ann* network),
                     double a, double b,double boundary_x,double boundary_y
        ,double boundary_ydot, double eps){
    int n = (network->params)->size;
    gsl_vector* step = gsl_vector_calloc(n);
    gsl_vector_memcpy(step,(network->params));
    for(int i=0; i<n; i++){
        double stepi = gsl_vector_get(step,i);
        gsl_vector_set(step,i,0.7*stepi);
    }
    N_MAX = 0;
    gsl_matrix* simplex = gsl_matrix_alloc(n,n+1);
    gsl_vector* centroid = gsl_vector_alloc(n);

    gsl_vector* init_vec = gsl_vector_alloc(n);
    ann* init = malloc(sizeof(ann));
    init->f = network->f;
    init->f_diff = network->f_diff;
    init->f_diffdiff = network->f_diffdiff;
    init->f_int = network->f_int;
    init->params = init_vec;

    ann* ann_p1 = malloc(sizeof(ann));
    gsl_vector* p1 = gsl_vector_alloc(n);
    ann_p1->f = (network->f);
    ann_p1->f_diff = (network->f_diff);
    ann_p1->f_diffdiff = (network->f_diffdiff);
    ann_p1->f_int = (network->f_int);
    ann_p1->params = p1;

    ann* ann_p2 = malloc(sizeof(ann));
    gsl_vector* p2 = gsl_vector_alloc(n);
    ann_p2->f = (network->f);
    ann_p2->f_diff = (network->f_diff);
    ann_p2->f_diffdiff = (network->f_diffdiff);
    ann_p2->f_int = (network->f_int);
    ann_p2 -> params = p2;

    gsl_vector* F_val = gsl_vector_alloc(n+1);
    int hi, lo;

    //starting point
    for(int i = 0; i<n+1; i++){
        gsl_matrix_set_col(simplex,i,(network->params));
    }
    for(int i = 0; i<n; i++){
        double ii = gsl_matrix_get(simplex,i,i);
        double stepi = gsl_vector_get(step,i);
        gsl_matrix_set(simplex,i,i, ii+stepi);
    }

    //finding high, low, centroid
    printf("init segmenttion?\n");
    annWild_init_vals(cost,simplex,init, diffeq_pow2, a, b, boundary_x, boundary_y
    ,boundary_ydot, F_val,centroid,&hi,&lo);
    //printf("%g\n",size(simplex,lo));
    //matrix_print(simplex,stdout);
    //printf("%g\n",size(simplex,lo));
    printf("while segmenttion?\n");
    while(size(simplex,lo)>eps && N_MAX<1e6){
        N_MAX++;
        hiLoCent(simplex,F_val,centroid,&hi,&lo);
        gsl_vector_view highest = gsl_matrix_column(simplex,hi);
        /*
        printf("N_MAX: %i \n ",N_MAX);

        printf("centroid: \n ");
        gsl_vector_fprintf(stdout,centroid,"%g");
        printf("Fvals: \n ");
        gsl_vector_fprintf(stdout,F_val,"%g");

        printf("Simplex: \n ");
        matrix_print(simplex,stdout);
        printf("highest = %i , lowest = %i\n", hi, lo);
        printf("Highest: \n");
        gsl_vector_fprintf(stdout,&highest.vector,"%g");
        */

        reflect(&highest.vector, centroid, p1);


        double f_re = cost(ann_p1, diffeq_pow2,
                           a, b, boundary_x, boundary_y
                ,boundary_ydot);

        //printf("P1: \n");
        //gsl_vector_fprintf(stdout, p1 ,"%g");

        if(f_re< gsl_vector_get(F_val,lo)){
            expand(&highest.vector,centroid,p2);
            double f_ex = cost(ann_p2,diffeq_pow2,
                               a, b, boundary_x, boundary_y
                    ,boundary_ydot);
            if(f_ex<f_re){
                gsl_vector_memcpy(&highest.vector,p2);
                gsl_vector_set(F_val,hi,f_ex);
            }
            else{
                gsl_vector_memcpy(&highest.vector,p1);
                gsl_vector_set(F_val,hi,f_re);
            }
        }
        else{
            if(f_re<gsl_vector_get(F_val,hi)){
                gsl_vector_memcpy(&highest.vector,p1);
                gsl_vector_set(F_val,hi,f_re);
            }

            else {
                contract(&highest.vector, centroid, p1);
                double f_co = cost(ann_p1,diffeq_pow2,
                                   a, b, boundary_x, boundary_y
                        ,boundary_ydot);
                if (f_co < gsl_vector_get(F_val, hi)) {
                    gsl_vector_memcpy(&highest.vector, p1);
                    gsl_vector_set(F_val, hi, f_co);
                }
                else {
                    reduce(simplex, lo);
                    annWild_init_vals(cost,simplex,init, diffeq_pow2,
                                  a, b, boundary_x, boundary_y
                            ,boundary_ydot ,F_val,centroid,&hi,&lo);
                }
            }

        }
    }
    gsl_vector_view final_low = gsl_matrix_column(simplex,lo);
    //printf("final\n");
    //gsl_vector_fprintf(stdout,&final_low.vector,"%g");
    gsl_vector_memcpy((network->params),&final_low.vector);
    //printf("params\n");
    //gsl_vector_fprintf(stdout,(network->params),"%g");
    printf("N_MAX = %i\n",N_MAX);
    gsl_vector_free(centroid);
    gsl_vector_free(F_val);
    gsl_vector_free(p1);
    gsl_vector_free(p2);
    free(ann_p1);
    free(ann_p2);
    free(init);
    gsl_matrix_free(simplex);
}
























//IF NESTED FUNCTIONS ARE ALLOWED ONE CAN USE THESE

void init_vals(double f(gsl_vector* x), gsl_matrix* simplex, gsl_vector* F_val, gsl_vector* centroid, int* hi, int* lo){
    int n = centroid -> size;
    for(int i = 0; i<n+1; i++){
        gsl_vector_view x = gsl_matrix_column(simplex,i);
        gsl_vector_set(F_val,i,f(&x.vector));
    }
    hiLoCent(simplex,F_val,centroid,hi,lo);
}

void amoeba( double f(gsl_vector* x),
             gsl_vector* x, gsl_vector* step, double eps){
    int n = x->size;
    N_MAX = 0;
    gsl_matrix* simplex = gsl_matrix_alloc(n,n+1);
    gsl_vector* centroid = gsl_vector_alloc(n);
    gsl_vector* p1 = gsl_vector_alloc(n);
    gsl_vector* p2 = gsl_vector_alloc(n);
    gsl_vector* F_val = gsl_vector_alloc(n+1);
    int hi, lo;

    //starting point
    for(int i = 0; i<n+1; i++){
        gsl_matrix_set_col(simplex,i,x);
    }
    for(int i = 0; i<n; i++){
        double ii = gsl_matrix_get(simplex,i,i);
        double stepi = gsl_vector_get(step,i);
        gsl_matrix_set(simplex,i,i, ii+stepi);
    }

    //finding high, low, centroid
    init_vals(f,simplex,F_val,centroid,&hi,&lo);
    while(size(simplex,lo)>eps && N_MAX<10000){
        N_MAX++;
        hiLoCent(simplex,F_val,centroid,&hi,&lo);
        gsl_vector_view highest = gsl_matrix_column(simplex,hi);
        /*
        printf("N_MAX: %i \n ",N_MAX);

        printf("centroid: \n ");
        gsl_vector_fprintf(stdout,centroid,"%g");

        printf("Simplex: \n ");
        matrix_print(simplex,stdout);
        printf("highest = %i , lowest = %i\n", hi, lo);
        printf("Highest: \n");
        gsl_vector_fprintf(stdout,&highest.vector,"%g");
        */
        reflect(&highest.vector, centroid, p1);


        double f_re = f(p1);
        //printf("P1: \n");
        //gsl_vector_fprintf(stdout, p1 ,"%g");

        if(f_re< gsl_vector_get(F_val,lo)){
            expand(&highest.vector,centroid,p2);
            double f_ex = f(p2);
            if(f_ex<f_re){
                gsl_vector_memcpy(&highest.vector,p2);
                gsl_vector_set(F_val,hi,f_ex);
            }
            else{
                gsl_vector_memcpy(&highest.vector,p1);
                gsl_vector_set(F_val,hi,f_re);
            }
        }
        else{
            if(f_re<gsl_vector_get(F_val,hi)){
                gsl_vector_memcpy(&highest.vector,p1);
                gsl_vector_set(F_val,hi,f_re);
            }

            else {
                contract(&highest.vector, centroid, p1);
                double f_co = f(p1);
                if (f_co < gsl_vector_get(F_val, hi)) {
                    gsl_vector_memcpy(&highest.vector, p1);
                    gsl_vector_set(F_val, hi, f_co);
                }
                else {
                    reduce(simplex, lo);
                    init_vals(f, simplex, F_val, centroid, &hi, &lo);
                }
            }

        }
    }
    gsl_vector_view final_low = gsl_matrix_column(simplex,lo);
    gsl_vector_memcpy(x,&final_low.vector);
    gsl_vector_free(centroid);
    gsl_vector_free(F_val);
    gsl_vector_free(p1);
    gsl_vector_free(p2);
    gsl_matrix_free(simplex);
}
