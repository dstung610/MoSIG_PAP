#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <math.h>

#define ndims 2
#define s 9
// Normal Matrix Multiplication
int main(int argc, char *argv[]){
        int size = s;
        if (argc >= 2) {
                size = atoi(argv[1]);
        }
        int matA[size], matB[size], matC[size];
        int dim  = sqrt(size);
        //Initializtion
        for (int i = 0; i < size; i++) {
                matA[i] = i;
                matB[i] = i;
                matC[i] = 0;
        }

        //Computation
        for (int j = 0; j < dim; j++) {
                for (int k = 0; k < dim; k++) {
                        for (int l = 0; l < dim; l++) {
                                matC[j*dim + k] += matA[j*dim + l] * matB[k+(l*dim)];
                        }
                }
        }

        //Display Result
        for (int m = 0; m < size; m++) {
                if (m%dim==0 && m!=0) {
                        printf("\n%d ", matC[m]);
                }
                else{
                        printf("%d ", matC[m]);
                }
        }
        printf("\n");
        return 0;
}
