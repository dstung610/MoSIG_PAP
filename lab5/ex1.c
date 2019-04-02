#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <math.h>

#define ndims 2

// Fox Algorithm Matrix Multiplication
int main(int argc, char *argv[]){
        int rank,size;
        int reorder = 0;
        int dims[ndims] = {0}, periods[ndims] = {0}, coords[ndims] = {0}, maxdims = ndims, remains[ndims]={0,1};
        int nnodes, dim;
        int a,b,c;
        MPI_Comm topologyComm, subComm; //subComm for broadcast A

        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        //Init value of matrix
        srand(time(NULL)+rank);
        // a = (rand()%10);
        // b = (rand()%10);
        a = rank;
        b = rank;

        dim = floor(sqrt(size));
        nnodes = dim*dim;

        for (int i = 0; i < ndims; i++) {
                dims[i] = dim;
        }

        printf("Rank %d ", rank);
        MPI_Dims_create(nnodes, ndims, dims);
        MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &topologyComm);
        if (rank < nnodes) {
                MPI_Cart_sub(topologyComm, remains, &subComm);
                MPI_Cart_coords(topologyComm, rank, maxdims, coords);
                printf("has coordinates (%d,%d) holds value A=%d, B=%d in topologyComm\n", coords[0], coords[1],a,b);
                MPI_Barrier(topologyComm);


                int rank_old = rank;
                int a_old, b_old, root;
                MPI_Comm_rank(subComm, &rank);
                MPI_Comm_size(subComm, &size);
                for (int stage = 0; stage < 1; stage++) {
                        a_old = a;
                        b_old = b;
                        if (coords[0]==coords[1]) {
                                root = ((rank_old%dim)+stage)%dim;
                        } else {

                        }

                        MPI_Bcast(&a,1,MPI_INT,root,subComm);
                        MPI_Barrier(subComm);
                        printf("\tStage %d Rank %d Root %d :A=%d B=%d => A=%d B=%d\n", stage, rank_old,((rank_old%dim)+stage)%dim, a_old, b_old, a, b);
                        MPI_Barrier(topologyComm);
                        a = a_old;
                }


        } else {
                printf("\n");
        }


        MPI_Finalize();
}
