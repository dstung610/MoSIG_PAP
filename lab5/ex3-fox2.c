#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <math.h>

#define ndims 2

void printMat(int mat[], int size, int dim){
        for (int i = 0; i < size; i++) {
                if (i%dim==0 && i!=0) {
                        printf("\n%d ", *(mat+i));
                } else {
                        printf("%d ", *(mat+i));
                }
        }
        printf("\n---------\n");
}

// Fox Algorithm Matrix Multiplication with data distribution
int main(int argc, char *argv[]){
        int rank_topology, rank_col_comm, rank_row_comm, size;
        int reorder = 0;
        int dims[ndims] = {0}, periods[ndims] = {1,1}, coords[ndims] = {0}, maxdims = ndims, remainsRow[ndims]={0,1}, remainsCol[ndims]={1,0};
        int nnodes, dim;
        int a = 0, b = 0, c = 0;
        MPI_Comm topologyComm, rowComm, colComm; //subComm for broadcast A and shifting B
        MPI_Status status;
        MPI_Datatype typeMatRow, typeMatCol;

        MPI_Init(&argc,&argv);

        MPI_Comm_rank(MPI_COMM_WORLD, &rank_topology);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        //Init value of matrix at rank 0
        dim = floor(sqrt(size));
        nnodes = dim*dim;

        MPI_Type_vector(dim, 1, 1, MPI_INT, &typeMatRow);
        MPI_Type_commit(&typeMatRow);
        MPI_Type_vector(dim, 1, dim, MPI_INT, &typeMatCol);
        MPI_Type_commit(&typeMatCol);

        int matA[nnodes], matB[nnodes], matC[nnodes], sendCount[nnodes], displacesRow[nnodes], displacesCol[nnodes];
        if (rank_topology == 0 ) {
                for (int i = 0; i < nnodes; i++) {
                        matA[i] = i;
                        matB[i] = i;
                        matC[i] = i;
                        sendCount[i] = 1;
                }
        } else {
                for (int i = 0; i < nnodes; i++) {
                        matA[i] = 0;
                        matB[i] = 0;
                        matC[i] = 0;
                        sendCount[i] = 0;
                }
        }


        for (int i = 0; i < nnodes; i++) {
                displacesRow[i] = i/dim;
                displacesCol[i] = i%dim;
        }

        // printMat(displacesCol,nnodes,dim);
        // printMat(displacesRow,nnodes,dim);

        for (int i = 0; i < ndims; i++) {
                dims[i] = dim;
        }

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Dims_create(nnodes, ndims, dims);
        MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &topologyComm);

        if (rank_topology < nnodes) {
                MPI_Cart_sub(topologyComm, remainsRow, &rowComm);
                MPI_Cart_sub(topologyComm, remainsCol, &colComm);
                MPI_Cart_coords(topologyComm, rank_topology, maxdims, coords);

                int a_temp, b_temp, offsetA, offsetB;
                // Data distribution
                // MPI_Scatterv(void *sendbuf,int *sendcnts,int *displs,MPI_Datatype sendtype,void *recvbuf,int recvcnt,MPI_Datatype recvtype,int root,MPI_Comm comm);
                MPI_Scatterv(matA, sendCount, displacesRow, typeMatRow, matA, 1, typeMatRow, 0, topologyComm);
                MPI_Scatterv(matB, sendCount, displacesCol, typeMatCol, matB, 1, typeMatRow, 0, topologyComm);

                printf("Rank %d: MatA %d %d %d - MatB %d %d %d ", rank_topology, matA[0],matA[1],matA[2],matB[0],matB[1],matB[2]);
                sleep(1);
                // Loop for n stages
                for (int stage = 0; stage < dim; stage++) {
                        offsetA = (coords[0]+stage)%dim;
                        offsetB = stage;
                        a_temp = matA[offsetA];
                        b_temp = matB[offsetB];

                        // Local Computation
                        c = c + (a_temp * b_temp);
                }
                printf("Rank %d has c = %d\n", rank_topology,c);
        } else {
                printf("\n");
        }


        MPI_Finalize();
        return 0;
}
