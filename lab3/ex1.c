#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]){
  int rank,size;
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("Rank:%d - Size:%d\n",rank, size);
  if (rank%2 == 0) {
    printf("I have rank %d and it's even\n", rank);
  }
  MPI_Finalize();
}
