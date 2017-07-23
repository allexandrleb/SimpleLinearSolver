#include <stdio.h>
#include <math.h>
#include "nlf.h"
#include "linalg.h"

#define SAMPLE_SIZE 25
#define PARAMS_DIM  7
#define NUM_PARAMS  7
#define FVAL_CONVERGE 1e-7
#define ITER_LIM 50

double eval_func(double x, const double *params){
    return params[0]+params[1]*tanh(params[2]*x+params[3])+params[4]*tanh(params[5]*(-x)+params[6]);
}

//derivative by x
double eval_deriv(double x, const double *params){
    //d(tanh(x))/dx = 1-(tanh(x))^2
    return  params[1]*params[2]*(1-pow(tanh(params[2]*x),2)) - params[3]*params[4]*(1-pow(tanh(params[4]*(-x)),2));
}

//gradient in space of parameters
void eval_grad(double x, const double *params,double *grad){

    //d(tanh(x))/dx = 1-(tanh(x))^2  = 1/(csh(x))^2

    grad[0] = 1.0;
    grad[1] = tanh(params[2]*x+params[3]);
    grad[2] = params[1]*(1-pow(tanh(params[2]*x+params[3]),2))*x;
    grad[3] = params[1]*(1-pow(tanh(params[2]*x+params[3]),2));
    grad[4] = tanh(params[5]*(-x)+params[6]);
    grad[5] = params[4]*(1-pow(tanh(params[5]*(-x)+params[6]),2))*(-x);
    grad[6] = params[4]*(1-pow(tanh(params[5]*(-x)+params[6]),2));

//    printf("gradient of f at x=%lf:\n",x);
//    print_colvec(grad);


}

//gradient of discrepancy function [F = \sum_i (y_i -f(x_i,p))]
void gradF(const double *x, const double *y, const double *params,double *Grad){

    //initialize with zeros
    for(size_t j=0; j < NUM_PARAMS; j++)
        Grad[j] = 0.0;

    for(size_t i=0; i < SAMPLE_SIZE; i++){
        double fgrad[NUM_PARAMS];
        eval_grad(x[i],params,fgrad);
    for(size_t j=0; j < NUM_PARAMS; j++)
        Grad[j] += -2.0* (y[i] - eval_func(x[i],params))*fgrad[j];
    }

}


//hessian of interpolation function
void eval_hessian(double x, const double *params,double *hessian){

    //second derivative with params[0]
    for(size_t i=0; i<PARAMS_DIM ; i++){
        hessian[i] = 0.0;
        hessian[i*PARAMS_DIM] = 0.0;
    }

    //second derivate with params[1]
    hessian[1+PARAMS_DIM*1] = 0.0;
    hessian[2+PARAMS_DIM*1] = (1.0-pow(tanh(params[2]*x+params[3]),2.0))*x;
    hessian[3+PARAMS_DIM*1] = (1.0-pow(tanh(params[2]*x+params[3]),2.0));
    hessian[4+PARAMS_DIM*1] = 0.0;
    hessian[5+PARAMS_DIM*1] = 0.0;
    hessian[6+PARAMS_DIM*1] = 0.0;
    for(size_t i=1; i<PARAMS_DIM ; i++)
        hessian[i*PARAMS_DIM+1] = hessian[i+PARAMS_DIM*1];


    //second derivate with params[2]
    hessian[2+PARAMS_DIM*2] = -2.0*params[1]*tanh(params[2]*x+params[3])*(1-pow(tanh(params[2]*x+params[3]),2.0))*x*x;
    hessian[3+PARAMS_DIM*2] = -2.0*params[1]*tanh(params[2]*x+params[3])*(1-pow(tanh(params[2]*x+params[3]),2.0))*x;
    hessian[4+PARAMS_DIM*2] = 0.0;
    hessian[5+PARAMS_DIM*2] = 0.0;
    hessian[6+PARAMS_DIM*2] = 0.0;
    for(size_t i=2; i<PARAMS_DIM ; i++)
        hessian[i*PARAMS_DIM+2] = hessian[i+PARAMS_DIM*2];


    //second derivate with params[3]
    hessian[3+PARAMS_DIM*3] = -2.0*params[1]*tanh(params[2]*x+params[3])*(1-pow(tanh(params[2]*x+params[3]),2.0));
    hessian[4+PARAMS_DIM*3] = 0.0;
    hessian[5+PARAMS_DIM*3] = 0.0;
    hessian[6+PARAMS_DIM*3] = 0.0;
    for(size_t i=3; i<PARAMS_DIM ; i++)
        hessian[i*PARAMS_DIM+3] = hessian[i+PARAMS_DIM*3];

    //second derivate with params[4]
    hessian[4+PARAMS_DIM*4] = 0.0;
    hessian[5+PARAMS_DIM*4] = (1-pow(tanh(params[5]*(-x)+params[6]),2.0))*(-x);
    hessian[6+PARAMS_DIM*4] = (1-pow(tanh(params[5]*(-x)+params[6]),2.0));
    for(size_t i=4; i<PARAMS_DIM ; i++)
        hessian[i*PARAMS_DIM+4] = hessian[i+PARAMS_DIM*4];


    //second derivate with params[5]
    hessian[5+PARAMS_DIM*5] = -2.0*params[4]*tanh(params[5]*(-x)+params[6])*(1-pow(tanh(params[5]*(-x)+params[6]),2.0))*(-x)*(-x);
    hessian[6+PARAMS_DIM*5] = -2.0*params[4]*tanh(params[5]*(-x)+params[6])*(1-pow(tanh(params[5]*(-x)+params[6]),2.0))*(-x);
    for(size_t i=5; i<PARAMS_DIM ; i++)
        hessian[i*PARAMS_DIM+5] = hessian[i+PARAMS_DIM*5];

    //second derivate with params[6]
    hessian[6+PARAMS_DIM*6] = -2.0*params[4]*tanh(params[5]*(-x)+params[6])*(1-pow(tanh(params[5]*(-x)+params[6]),2.0));

//    printf("hessian of f at x=%lf and params: ",x);
//    print_rowvec(params);
//    print_matr(hessian);



}

//hessian of discrepancy function [F = \sum_i (y_i -f(x_i,p))]
void hessF(const double *x, const double *y, const double *params,double *Hess){

    //initialize with zeros
    for(size_t j=0; j < NUM_PARAMS*NUM_PARAMS; j++)
        Hess[j] = 0.0;

    for(size_t i=0; i < SAMPLE_SIZE; i++){

    double fgrad[NUM_PARAMS];
    eval_grad(x[i],params,fgrad);

    double fhess[NUM_PARAMS*NUM_PARAMS];
    eval_hessian(x[i], params, fhess);

    //Hess is actually symmetric, this can be used to save computation
    for(size_t j=0; j < NUM_PARAMS; j++){
    for(size_t k=0; k < NUM_PARAMS; k++){
        Hess[j+NUM_PARAMS*k] += 2.0*fgrad[j]*fgrad[k] - 2.0 * (y[i] - eval_func(x[i],params))*fhess[j+NUM_PARAMS*k] ;
    }
    }
    }


} 






void gen_fit(const double *x,const double *y){

    fit_conjugate_gradient(x,y);

    fit_steepest_descent(x,y);
    fit_gauss_newton(x,y);
    fit_levenberg_marquardt(x,y);
    fit_linear_system(x,y);

    printf("fitting performed\n");

}


void fit_conjugate_gradient(const double *x,const double *y){

    printf("fitting by conjugate gradient method complete\n");

}

void fit_steepest_descent(const double *x,const double *y){

    //double params[PARAMS_DIM]={2.50, 2.10, 2.00, 0.00, -1.40, 1.00, 0.00};
    double params[PARAMS_DIM]={3.53, -2.12, 4.03, 0.00, -1.48, 1.00, 0.00};
    int optparams[PARAMS_DIM]={1, 1, 1, 0, 1, 0, 0};
    double steps[PARAMS_DIM] ={0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00};
    double residuals[SAMPLE_SIZE];
    size_t count=0;
    double F=1.0;
    double Fprev=1.0;
    double gradF[PARAMS_DIM]={0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};



    while(F>FVAL_CONVERGE || count<1){
        Fprev=F;
        F=0.0;
#pragma omp parallel for
        for(size_t j=0; j<NUM_PARAMS; j++)
            gradF[j]=0.0;

#pragma omp parallel for
        for(size_t i=0;i<SAMPLE_SIZE; i++){
            residuals[i] = eval_func(x[i],params)-y[i];
            F+=residuals[i]*residuals[i];
            double g[PARAMS_DIM];
            eval_grad(x[i],params,g);
            for(size_t j=0; j<NUM_PARAMS; j++){
                if(optparams[j])
                    gradF[j]+= 2.0*residuals[i]*g[j];
            }
        }
        double gradNorm  = 0.0;
 //       double pstepNorm  = 0.0;
 //       double gradProj   = 0.0;
 //       double maxGrad = 0.0;
#pragma omp parallel for
        for(size_t j=0; j<NUM_PARAMS; j++){
            gradNorm +=gradF[j]*gradF[j];
//            pstepNorm +=steps[j]*steps[j];
//            gradProj +=steps[j]*gradF[j];
        }
        int stepdim=0;
//#pragma omp parallel for
//        for(size_t j=0; j<NUM_PARAMS; j++){
//            if(pstepNorm > 0){
//
//                gradF[j] -= gradProj/sqrt(pstepNorm)*steps[j];
//            }
//
//            if(fabs(gradF[j])>maxGrad) maxGrad=fabs(gradF[j]);
//        }

#pragma omp parallel for
        for(size_t j=0; j<NUM_PARAMS; j++)
           if(optparams[j]) stepdim++; 

        gradNorm = sqrt(gradNorm);
#pragma omp parallel for
        for(size_t j=0; j<NUM_PARAMS; j++)
//            steps[j] = (fabs(gradF[j])>0.1*maxGrad) ?  F/(gradF[j]+steps[j])/stepdim : 0.0;
            steps[j] =  (fabs(gradF[j]/(F+1e-15)) > 1.0) ? F/(gradF[j])/stepdim : 1.0/stepdim*gradF[j]/fabs(gradF[j]+1e-15) ;
//            steps[j] = (fabs(gradF[j])>1e-15) ? F/(gradF[j]) : 0.0;
        count++;

        //output
//        printf("iter %d: F = %g, |gradF| = %g\n",count,F,gradNorm);
//        printf("iter %d: %3s %10s %10s %10s\n",count,"[i]","param[i]","step[i]","gradF[i]");
//        for(size_t k=0; k< PARAMS_DIM; k++)
//            printf("iter %d: [%1d] %8.4g %10.3g %10.3g\n",count,k,params[k],steps[k],gradF[k]);
//
//        for(size_t k=0; k< 50; k++)
//            printf("-");
//        printf("\n");

        for(size_t j=0; j<NUM_PARAMS; j++) 
            params[j] -= steps[j];

        if(count > ITER_LIM){
            printf("OUT OF ITERATIONS LIMIT. ABORT\n");
            break;
        }
    }

    printf("fitting by steepest descent method complete\n");
    printf("optimized values of parameters:\n");
    for(size_t k=0; k< PARAMS_DIM; k++)
        printf("%6.3f ",params[k]);
    printf("\n");

}

void fit_gauss_newton(const double *x,const double *y){

    printf("fitting by gauss newton method complete\n");

}

void fit_levenberg_marquardt(const double *x,const double *y){

    printf("fitting by levenberg marquardt method complete\n");

}

void fit_linear_system(const double *x,const double *y){
//    double params[PARAMS_DIM]={3.53, -2.12, 4.03, 0.00, -1.48, 1.00, 0.00};
    double params[PARAMS_DIM]={3.53, 2.12, 4.03, 0.00, -1.48, 1.00, 0.00};
    double step[PARAMS_DIM]={1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
    int count = 0;

    while( (veclen(step) > 1e-5) && count < ITER_LIM){
    double FGrad[PARAMS_DIM];
    gradF(x,y,params,FGrad);


    for(size_t k=0; k< 50; k++)
            printf("-");
        printf("\n");



    printf("Residual:\n");
    print_colvec(FGrad);


    double FHess[PARAMS_DIM*PARAMS_DIM];
    hessF(x,y,params,FHess);


    printf("Hessian:\n");
    print_matr(FHess);


    double Jinv[PARAMS_DIM*PARAMS_DIM];
    invert(FHess,Jinv);


    printf("Inverse of Hessian:\n");
    print_matr(Jinv);

    for(size_t k=0; k< 50; k++)
        printf("-");
    printf("\n");



    mulcolvec(Jinv,FGrad,step);
    for(size_t i =0 ; i <PARAMS_DIM; i++)
        params[i] -= step[i];

    count++;
    }

    printf("params:\n");
    for(size_t i =0 ; i <PARAMS_DIM; i++)
        printf("%6.3f ",params[i]);
    printf("\n");

    printf("fitting by linear system is complete\n");

}
