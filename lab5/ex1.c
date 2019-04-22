#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <math.h>

#define ndims 2

// Fox Algorithm Matrix Multiplication
int main(int argc, char *argv[]){
        int rank_topology, rank_col_comm, rank_row_comm, size;
        int reorder = 0;
        int dims[ndims] = {0}, periods[ndims] = {1,1}, coords[ndims] = {0}, maxdims = ndims, remainsRow[ndims]={0,1}, remainsCol[ndims]={1,0};
        int nnodes, dim;
        int a, b, c = 0;
        MPI_Comm topologyComm, rowComm, colComm; //subComm for broadcast A and shifting B
        MPI_Status status;

        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_topology);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        //Init value of matrix
        srand(time(NULL)+rank_topology);
        // a = (rand()%10);
        // b = (rand()%10);
        a = rank_topology;
        b = rank_topology;

        dim = floor(sqrt(size));
        nnodes = dim*dim;

        for (int i = 0; i < ndims; i++) {
                dims[i] = dim;
        }


        MPI_Dims_create(nnodes, ndims, dims);
        MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &topologyComm);
        if (rank_topology < nnodes) {
                MPI_Cart_sub(topologyComm, remainsRow, &rowComm);
                MPI_Cart_sub(topologyComm, remainsCol, &colComm);
                MPI_Cart_coords(topologyComm, rank_topology, maxdims, coords);
                printf("Rank %d has coordinates (%d,%d) holds value A=%d, B=%d in topologyComm\n", rank_topology, coords[0], coords[1],a,b);

                int a_temp, root, src, dest, b_temp = b;
                MPI_Comm_rank(rowComm, &rank_row_comm);

                // Loop for n stages
                for (int stage = 0; stage < dim; stage++) {
                        a_temp = a;
                        b_temp = b;

                        // Broadcasting A
                        root = (coords[0] + stage) % dim;
                        MPI_Bcast(&a_temp,1,MPI_INT,root,rowComm);

                        // Local Computation
                        c = c + (a_temp * b_temp);

                        // Shifting B
                        MPI_Cart_shift(colComm, 0, -1, &src, &dest);
                        MPI_Sendrecv(&b_temp, 1, MPI_INT,dest, 95, &b, 1, MPI_INT, src, 95, colComm,&status);
                        MPI_Barrier(colComm);

                        //printf("\tStage %d Rank %d Root %d: A=%d B=%d\n", stage, rank_topology, root, a_temp, b_temp);
                }
                printf("Rank %d has c = %d\n", rank_topology,c);
        } else {
                printf("\n");
        }


        MPI_Finalize();
        return 0;
}
