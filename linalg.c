#include <stdio.h>
#include <math.h>
#include "linalg.h"

//#define MATR_WIDTH   7
//#define MATR_HEIGHT  7
#define MATR_N_ELEMS 49

#define MYZEROVAL 1e-10

#define DIVIDE(a,b) (a)/(MYZEROVAL+b)

void print_matr(const double* matr){
    for(size_t i = 0; i< MATR_HEIGHT; i++){
        for(size_t j = 0; j< MATR_WIDTH; j++){
            printf("%6.3lf ",matr[i*MATR_WIDTH+j]);
        }
        printf("\n");
    }
}

void print_colvec(const double* colvec){
    for(size_t i = 0; i< MATR_HEIGHT; i++){
        printf("%6.3lf\n",colvec[i]);
    }

}

void print_rowvec(const double* rowvec){
    for(size_t j = 0; j< MATR_WIDTH; j++){
        printf("%6.3lf ",rowvec[j]);
    }
    printf("\n");

}

double veclen(const double* vec){

    double len = 0.0;
    for(size_t i =0 ; i <MATR_HEIGHT; i++)
        len += vec[i]*vec[i];
    return sqrt(len);

}

void multiply(const double* L,const double* R, double* Out){

    for(size_t i=0;i<MATR_HEIGHT;i++){
    for(size_t j=0;j<MATR_WIDTH;j++){
        size_t idxo = i*MATR_WIDTH+j;
        Out[idxo] = 0.0;
    for(size_t k=0;k<MATR_WIDTH;k++){
        size_t idxl = i*MATR_WIDTH+k;
        size_t idxr = k*MATR_WIDTH+j;
        Out[idxo] += L[idxl]*R[idxr];
    }
    }
    }
}



void gen_reflection(int kpiv, const double* V, double* H){

    //V - vector used to generate reflection
    //kpiv - index of pivot element of vector
    //H - output parameter, the reflection matrix

    //double V[MATR_HEIGHT] = { 1.0, 2.0, -3.0, -2.0, 3.0,0.0,0.0};
    //int kpiv = 2;
    int sign = V[kpiv] > 0  ? +1 : -1;

    
    double h[MATR_HEIGHT]; //difference vector between the source and the target
    for(size_t i =0 ; i <MATR_HEIGHT; i++)
        h[i] = V[i];

    double lenV=0.0;
    for(size_t i=kpiv;i<MATR_HEIGHT; i++)
        lenV += V[i]*V[i];
    lenV=sqrt(lenV);

    double ek[MATR_HEIGHT];
    for(size_t i =0 ; i <MATR_HEIGHT; i++){
        ek[i] = (i==kpiv) ? 1.0 : 0.0;
    }
    
    for(size_t i =kpiv ; i <MATR_HEIGHT; i++){
         h[i] -= (double)sign*lenV * ek[i];
    }

    double lenh=0.0;
    for(size_t i=kpiv;i<MATR_HEIGHT; i++)
        lenh += h[i]*h[i];
    lenh=sqrt(lenh);
    if(lenh>1e-15){
        for(size_t i=kpiv;i<MATR_HEIGHT; i++)
            h[i] = h[i]/lenh;
    }



//    double H[MATR_N_ELEMS]; //reflection matrix
    for(size_t i=0; i<MATR_HEIGHT; i++){
    for(size_t j=0; j<MATR_WIDTH; j++){
        size_t idx = i*MATR_WIDTH + j;
        if( i==j) H[idx]=1.0;
        else      H[idx]=0.0;
        if(i>=kpiv && j>=kpiv)
            H[idx] -= 2*h[i]*h[j];
    }
    }

    //printf("column:\n");
    //print_colvec(V);
    //printf("resulting reflection matrix for the column\n");
    //print_matr(H);


}
void transpose(const double* H, double* Q){
    for(size_t j =0; j< MATR_WIDTH; j++){
    for(size_t i = 0; i < MATR_HEIGHT; i++)
        Q[i*MATR_WIDTH+j] = H[j*MATR_WIDTH+i];
    }
}

void gen_reflection_by_matr(const double* M, double* H){
    //M - Matrix of columns to generate transform
    //H - orthogonal transform to triangulaze matrix M

    //double H[MATR_N_ELEMS];
    double Mtmp[MATR_N_ELEMS];
    for(size_t j =0; j< MATR_WIDTH; j++){
    for(size_t i = 0; i < MATR_HEIGHT; i++)
        Mtmp[i*MATR_WIDTH+j] = M[i*MATR_WIDTH+j];
    }

    //initial H is identity matrix
    for(size_t j =0; j< MATR_WIDTH; j++){
    for(size_t i = 0; i < MATR_HEIGHT; i++)
        H[i*MATR_WIDTH+j] = i==j ? 1.0 : 0.0;
    }

    for(size_t j =0; j< MATR_WIDTH; j++){
    double Col[MATR_HEIGHT];
    for(size_t i = 0; i < MATR_HEIGHT; i++)
        Col[i] = Mtmp[i*MATR_WIDTH+j];
    double Hj[MATR_N_ELEMS];
    gen_reflection(j,Col,Hj);
    double Tmp[MATR_N_ELEMS];
    multiply(Hj,H,Tmp);
        //overwrite H
    for(size_t j =0; j< MATR_WIDTH; j++){
    for(size_t i = 0; i < MATR_HEIGHT; i++)
        H[i*MATR_WIDTH+j] = Tmp[i*MATR_WIDTH+j];
    }

    multiply(Hj,Mtmp,Tmp);
        //overwrite Mtmp
    for(size_t j =0; j< MATR_WIDTH; j++){
    for(size_t i = 0; i < MATR_HEIGHT; i++)
        Mtmp[i*MATR_WIDTH+j] = Tmp[i*MATR_WIDTH+j];
    }



    }
}

void invert_triangle(const double *R, double *Inv){

    double Rtmp[MATR_N_ELEMS];
    //initialize Inv as identity matrix
    for(size_t j =0; j< MATR_WIDTH; j++){
    for(size_t i = 0; i < MATR_HEIGHT; i++){
        Inv[i*MATR_WIDTH+j] = i==j ? 1.0 : 0.0;
        Rtmp[i*MATR_WIDTH+j] = R[i*MATR_WIDTH+j];
    }
    }

    //
    //gaussian elimination:
    //Inv(ij) = ( I(ij) - \sum_{k>i} R(i,k)*Inv(k,j) ) / R(i,i)
    //
    for(int i = MATR_HEIGHT-1 ; i >=0; i--){
    for(int j = i ; j < MATR_WIDTH; j++){

        //scaling of original values at position (i,j)
        Inv[i*MATR_WIDTH+j] = DIVIDE(Inv[i*MATR_WIDTH+j],R[i*MATR_WIDTH+i]);
        Rtmp[i*MATR_WIDTH+j] = DIVIDE(Rtmp[i*MATR_WIDTH+j],R[i*MATR_WIDTH+i]);

        //subtraction of lower rows 
        for(int k = MATR_HEIGHT-1 ; k >i; k--){
            Inv[i*MATR_WIDTH+j] -= DIVIDE(R[i*MATR_WIDTH+k],R[i*MATR_WIDTH+i]) * Inv[k*MATR_WIDTH+j];
            Rtmp[i*MATR_WIDTH+j] -= DIVIDE(R[i*MATR_WIDTH+k],R[i*MATR_WIDTH+i])* Rtmp[k*MATR_WIDTH+j];
        }
    }

    //print progress of gaussian elimination
    //for(size_t ii = 0; ii < MATR_HEIGHT; ii++){
    //    for(size_t jj =0; jj< MATR_WIDTH; jj++)
    //        printf("%5.2f ",Rtmp[ii*MATR_WIDTH+jj]);
    //    printf(" | ");
    //    for(size_t jj =0; jj< MATR_WIDTH; jj++)
    //        printf("%5.2f ",Inv[ii*MATR_WIDTH+jj]);
    //    printf("\n");

    //}
    //for(size_t jj =0; jj< MATR_WIDTH*12+3; jj++)
    //    printf("-");
    //printf("\n");

    }


}

void mulcolvec(const double *M, const double *V, double *Res){

    for(size_t i = 0; i < MATR_HEIGHT; i++)
        Res[i] = 0.0;

    for(size_t i = 0; i < MATR_HEIGHT; i++){
    for(size_t j =0; j< MATR_WIDTH; j++)
        Res[i] = M[i*MATR_WIDTH+j]*V[j];
    }

}


void invert(const double *M,const double *Minv){

    double H[MATR_N_ELEMS];
    gen_reflection_by_matr(M,H);
    double R[MATR_N_ELEMS];

    multiply(H,M,R);

    double Q[MATR_N_ELEMS];
    transpose(H,Q);

    double Rinv[MATR_N_ELEMS];
    invert_triangle(R,Rinv);

    double Qinv[MATR_N_ELEMS];
    transpose(Q,Qinv);

   // double Minv[MATR_N_ELEMS];
    multiply(Rinv,Qinv,Minv);

}


//int main(){
int do_linalg(){

    double M[MATR_N_ELEMS]= { 1.0, 2.0, 0.0, 3.0, 5.0, 1.0 ,1.0,
                              0.0, 1.0,-1.0, 2.0, 4.0, 1.0 ,1.0,
                             -2.0, 0.0, 2.0, 4.0, 0.0, 1.0 ,1.0,
                              0.0,-3.0, 1.0, 3.0, 0.0, 1.0 ,1.0,
                              0.0, 0.0, 4.0, 1.0, 2.0, 3.0 ,1.0,
                              0.0, 0.0, 1.0, 1.0, 2.0, 1.0 ,1.0,
                              0.0, 2.0, 0.0, 1.0, 2.0,-1.0 ,1.0};

    double H[MATR_N_ELEMS];
    gen_reflection_by_matr(M,H);
    double R[MATR_N_ELEMS];

    multiply(H,M,R);

    double Q[MATR_N_ELEMS];
    transpose(H,Q);

    double Mr[MATR_N_ELEMS];
    multiply(Q,R,Mr);

    printf("orthogonal matrix:\n");
    print_matr(Q);
    printf("triangular matrix:\n");
    print_matr(R);
    printf("reconstructed matrix:\n");
    print_matr(Mr);
    printf("original matrix:\n");
    print_matr(M);

    double Rinv[MATR_N_ELEMS];
    invert_triangle(R,Rinv);

    //printf("inverted triangular matrix:\n");
    //print_matr(Rinv);


    //double Rr[MATR_N_ELEMS];
    //multiply(Rinv,R,Rr);
    //printf("invert validation:\n");
    //print_matr(Rr);

    double Qinv[MATR_N_ELEMS];
    transpose(Q,Qinv);

    double Minv[MATR_N_ELEMS];
    multiply(Rinv,Qinv,Minv);

    printf("inverted Full matrix M:\n");
    print_matr(Minv);

    multiply(Minv,M,Mr);
    printf("invert validation:\n");
    print_matr(Mr);





//    double Row[MATR_WIDTH]  = {1.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0};
//
//    double Col[MATR_HEIGHT] = {1.0, 1.0, 0.0, 0.0, 0.0, 2.0, 2.0};
//
//    printf("Matrix M:\n");
//    print_matr(M);
//
//    printf("Row-vector Row:\n");
//    print_rowvec(Row);
//
//    printf("Col-vector Col:\n");
//    print_colvec(Col);

//    double b[MATR_HEIGHT] = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//    double h[MATR_HEIGHT];
//    double lensq=0.0;
//    double lenCsq=0.0;
//    for(size_t i=0; i< MATR_HEIGHT; i++)
//        lenCsq += C[i] * C[i];
//
//    for(size_t i=0; i< MATR_HEIGHT; i++){
//        h[i] =  (sqrt(lenCsq)*b[i] - C[i]);
//        lensq += h[i] * h[i];
//    }
//    for(size_t i=0; i< MATR_HEIGHT; i++){
//        h[i] = h[i] /sqrt(lensq);
//    }
//    printf("Col-vector h:\n");
//    print_colvec(h);

//    double H[MATR_N_ELEMS];
//    for(size_t i=0; i<MATR_HEIGHT; i++){
//    for(size_t j=0; j<MATR_WIDTH; j++){
//        size_t idx = i*MATR_WIDTH + j;
//        if( i==j) H[idx]=1.0;
//        else      H[idx]=0.0;
//        H[idx] -= 2*h[i]*h[j];
//    }
//    }
//
//    printf("Matrix H:\n");
//    print_matr(H);
//
//    double Cr[MATR_HEIGHT];
//    for(size_t i =0; i< MATR_HEIGHT; i++){
//        Cr[i] = 0.0;
//    for(size_t j =0; j< MATR_WIDTH; j++){
//        Cr[i] += H[i*MATR_WIDTH+j]*C[j];
//    }
//    }
//
//    printf("Reflected C:\n");
//    print_colvec(Cr);
//
//    Numerical Analysis, Burden and Faires, 8th Edition

//    double V[MATR_HEIGHT] = { 1.0, 2.0, -3.0, -2.0, 3.0, 1.0, 2.0};
//    int kpiv = 2;
//    int sign = V[kpiv] > 0  ? +1 : -1;
//
//    
// //   double h[MATR_HEIGHT]; //difference vector between the source and the target
//    for(size_t i =0 ; i <MATR_HEIGHT; i++)
//        h[i] = V[i];
//
//    double lenV=0.0;
//    for(size_t i=kpiv;i<MATR_HEIGHT; i++)
//        lenV += V[i]*V[i];
//    lenV=sqrt(lenV);
//
//    double ek[MATR_HEIGHT];
//    for(size_t i =0 ; i <MATR_HEIGHT; i++){
//        ek[i] = (i==kpiv) ? 1.0 : 0.0;
//    }
//    
//    for(size_t i =kpiv ; i <MATR_HEIGHT; i++){
//         h[i] -= (double)sign*lenV * ek[i];
//    }
//
//    double lenh=0.0;
//    for(size_t i=kpiv;i<MATR_HEIGHT; i++)
//        lenh += h[i]*h[i];
//    lenh=sqrt(lenh);
//    for(size_t i=kpiv;i<MATR_HEIGHT; i++)
//        h[i] = h[i]/lenh;
//
//
//
////    double H[MATR_N_ELEMS]; //reflection matrix
//    for(size_t i=0; i<MATR_HEIGHT; i++){
//    for(size_t j=0; j<MATR_WIDTH; j++){
//        size_t idx = i*MATR_WIDTH + j;
//        if( i==j) H[idx]=1.0;
//        else      H[idx]=0.0;
//        if(i>=kpiv && j>=kpiv)
//            H[idx] -= 2*h[i]*h[j];
//    }
//    }
//
//
//    double Vr[MATR_HEIGHT];
//    for(size_t i =0; i< MATR_HEIGHT; i++){
//        Vr[i] = 0.0;
//    for(size_t j =0; j< MATR_WIDTH; j++){
//        Vr[i] += H[i*MATR_WIDTH+j]*V[j];
//    }
//    }
//
//    printf("Source vector V:\n");
//    print_colvec(V);
//
//
//
//    printf("Matrix H:\n");
//    print_matr(H);
//
//
//    printf("Reflected vector Vr:\n");
//    print_colvec(Vr);
//
//

    return 0;
}
