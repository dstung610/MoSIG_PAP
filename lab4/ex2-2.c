#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]){
  int rank,size;
  int total;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Allreduce(&rank, &total, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  printf("Total is: %d\n", total);
  MPI_Finalize();
}
