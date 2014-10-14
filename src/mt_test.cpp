#include "FRandom.h"
#include <stdio.h>

int main(int arg_c, char** args){
    int N = 1000;
    double counts = 1e8;
    double* x = new double[N];
    double* y = new double[N];
    double dx = 1.0/N;
    FRandom ng(1);
    for(int i = 0; i<N; i++){
        x[i] = (i+1)*dx;
    }
    for(int z= 0; z<counts; z++){

        double v = ng.nextDouble();

        int dex = v==1.0?N-1:int(v/dx);

        if( (dex<0) || (dex>=N) ){
            printf("failed index out of range %d\n", dex);
            throw "negative number";
        }

        y[dex] += 1;


    }

    for(int i = 0; i<N; i++){
        printf("%f\t%f\n", x[i], y[i]/counts);
    }
}