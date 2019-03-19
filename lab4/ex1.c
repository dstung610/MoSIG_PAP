#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]){
  int rank,size;
  int data[10];
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (rank==0) {
    data[1] = 5;
  }
  MPI_Bcast(data,10,MPI_INT,0,MPI_COMM_WORLD);
  if (rank == 0) {
    printf("I received data[1]: %d\n",data[1]);
  }
  MPI_Finalize();
}
