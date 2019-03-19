#include <stdio.h>
#include <mpi.h>
#include <unistd.h>

#define ndims 2

int main(int argc, char *argv[]){
  int rank,size;
  int dims[ndims];
  int nnodes;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // int MPI_Dims_create(int nnodes, int ndims, int dims[])
  // nnodes: number of nodes in a grid (integer)
  // ndims: number of cartesian dimensions (integer)
  // dims: integer array of size ndims specifying the number of nodes
  // in each dimension. A value of 0 indicates that MPI_Dims_create
  // should fill in a suitable value.
  MPI_Dims_create(nnodes, ndims, dims);

  // int MPI_Cart_create(MPI_Comm comm_old, int ndims, const int dims[],
  //                     const int periods[], int reorder, MPI_Comm * comm_cart)
  // Input Parameters
  // comm_old: input communicator (handle)
  // ndims: number of dimensions of cartesian grid (integer)
  // dims: integer array of size ndims specifying the number of processes in each dimension
  // periods: logical array of size ndims specifying whether the grid is periodic (true) or not (false) in each dimension
  // reorder: ranking may be reordered (true) or not (false) (logical)
  // Output Parameters
  // comm_cart: communicator with new cartesian topology (handle)
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims,);

  // int MPI_Cart_coords(MPI_Comm comm, int rank, int maxdims, int coords[])
  // Input Parameters
  // comm: communicator with cartesian structure (handle)
  // rank: rank of a process within group of comm (integer)
  // maxdims:length of vector coords in the calling program (integer)
  // Output Parameters
  // coords: integer array (of size ndims) containing the Cartesian coordinates of specified process (integer)
  MPI_Cart_coords();
  MPI_Finalize();
}
