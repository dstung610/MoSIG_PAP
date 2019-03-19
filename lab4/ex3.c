#include <stdio.h>
#include <mpi.h>
#include <unistd.h>

int main(int argc, char *argv[]){
  int rank,size;
  int data[10];
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  sleep(rank);
  printf("Rank %d is waiting...\n", rank);
  MPI_Barrier(MPI_COMM_WORLD);
  printf("Rank %d finished\n", rank);
  MPI_Finalize();
}
