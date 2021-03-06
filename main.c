#include <stdio.h>
#include <math.h>
#include "nlf.h"

int main(){


//    double x[25] = { 0.0,  1.0,  2.0,  3.0,  4.0,
//                     5.0,  6.0,  7.0,  8.0,  9.0,
//                    10.0, 11.0, 12.0, 13.0, 14.0,
//                    15.0, 16.0, 17.0, 18.0, 19.0,
//                    20.0, 21.0, 22.0, 23.0, 24.0};
//
//    double y[25] = { 0.0,  1.0,  2.0,  3.0,  4.0,
//                     5.0,  6.0,  7.0,  8.0,  9.0,
//                    10.0, 11.0, 12.0, 13.0, 14.0,
//                    15.0, 16.0, 17.0, 18.0, 19.0,
//                    20.0, 21.0, 22.0, 23.0, 24.0};

    double x[25];
    double y[25];

    for(size_t i =0 ; i <25; ++i){
        x[i] = (double)i - 5.0;
        y[i] = 2.5 + 2.1*tanh(2.0*x[i])-1.4*tanh(1.0*(-x[i]));
    }


    gen_fit(x,y);

    return 0;
}
