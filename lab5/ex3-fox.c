#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <math.h>

#define ndims 2

void printMat(int *mat, int size){
        int dim = sqrt(size);
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
        int *matA, *matB, *matC, *temp, *sendCount, recvCount, *displacesRow, *displacesCol;
        MPI_Comm topologyComm, rowComm, colComm; //subComm for broadcast A and shifting B
        MPI_Status status;
        MPI_Datatype typeMatRow, typeMatCol;



        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_topology);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        //Init value of matrix at rank 0
        dim = floor(sqrt(size));
        nnodes = dim*dim;

        if (rank_topology == 0) {
                matA = malloc(sizeof(int)*nnodes);
                matB = malloc(sizeof(int)*nnodes);
                matC = malloc(sizeof(int)*nnodes);

                for (int l = 0; l < nnodes; l++) {
                        temp = matA + l;
                        *temp = l;
                }

                for (int l = 0; l < nnodes; l++) {
                        temp = matB + l;
                        *temp = l;

                }
                printMat(matA,nnodes);
                printMat(matB,nnodes);
        }
        else{
                matA = malloc(sizeof(int)*dim);
                matB = malloc(sizeof(int)*dim);
                matC = malloc(sizeof(int)*dim);

                for (int l = 0; l < nnodes; l++) {
                        temp = matA + l;
                        *temp = 0;
                        temp = matB + l;
                        *temp = 0;
                }
        }

        recvCount = dim;
        sendCount = malloc(sizeof(int)*nnodes);
        displacesRow = malloc(sizeof(int)*nnodes);
        displacesCol = malloc(sizeof(int)*nnodes);

        for (int i = 0; i < nnodes; i++) {
                temp = sendCount + i;
                *temp = 1;
        }

        int offset = 0;
        for (int i = 0; i < nnodes; i++) {
                offset = i/dim;
                temp = displacesRow + i;
                *temp = offset;
        }

        offset = 0;
        for (int i = 0; i < nnodes; i++) {
                offset = i%dim;
                temp = displacesCol + i;
                *temp = offset;
        }

        if(rank_topology == 0 && 0 == 1) {
                printMat(sendCount,nnodes);
                printf("------\n");
                printMat(displacesRow,nnodes);
                printf("------\n");
                printMat(displacesCol,nnodes);
        }

        for (int i = 0; i < ndims; i++) {
                dims[i] = dim;
        }

        MPI_Type_vector(dim, 1, 1, MPI_INT, &typeMatRow);
        MPI_Type_vector(dim, 1, dim, MPI_INT, &typeMatCol);
        MPI_Type_commit(&typeMatRow);
        MPI_Type_commit(&typeMatCol);

        // if (rank_topology == 0) {
        //         MPI_Send(matB+2, 1, typeMatCol, 1, 95, MPI_COMM_WORLD);
        // }
        // if (rank_topology == 1) {
        //         MPI_Recv(matB, 1, typeMatCol, 0,95, MPI_COMM_WORLD, &status);
        //         printMat(matB,nnodes);
        // }
        //
        // sleep(2);
        // return 0;

        MPI_Dims_create(nnodes, ndims, dims);
        MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &topologyComm);
        MPI_Barrier(topologyComm);
        if (rank_topology < nnodes) {
                MPI_Cart_sub(topologyComm, remainsRow, &rowComm);
                MPI_Cart_sub(topologyComm, remainsCol, &colComm);
                MPI_Cart_coords(topologyComm, rank_topology, maxdims, coords);
                // printf("Rank %d has coordinates (%d,%d) holds value A=%d, B=%d in topologyComm\n", rank_topology, coords[0], coords[1],a,b);

                int a_temp, root, src, dest, b_temp, offsetA, offsetB;
                MPI_Comm_rank(rowComm, &rank_row_comm);

                // Data distribution
                // MPI_Scatterv(void *sendbuf,int *sendcnts,int *displs,MPI_Datatype sendtype,void *recvbuf,int recvcnt,MPI_Datatype recvtype,int root,MPI_Comm comm);
                MPI_Scatterv(matA, sendCount, displacesRow, typeMatRow, matA, 1, typeMatRow, 0, topologyComm);
                MPI_Barrier(topologyComm);
                MPI_Scatterv(matB, sendCount, displacesCol, typeMatCol, matB, 1, typeMatRow, 0, topologyComm);
                MPI_Barrier(topologyComm);
                printf("Rank %d MatA: %d %d %d\n", rank_topology,*matA, *(matA+1), *(matA+2));
                printf("Rank %d MatB: %d %d %d\n", rank_topology, *matB, *(matB+1), *(matB+2));


                // Loop for n stages
                for (int stage = 0; stage < dim; stage++) {
                        offsetA = (coords[0]+stage)%dim;
                        offsetB = stage;
                        a_temp = *(matA+offsetA);
                        b_temp = *(matB+offsetB);

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
