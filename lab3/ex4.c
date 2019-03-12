#include <stdio.h>

#include <mpi.h>

#define SIZE 10

void print_array(int array[]){
  for (int i = 0; i < SIZE; i++) {
  }
  printf("%2.f-", array[i]);
  printf("\n");
}

int main(int argc, char *argv[]){
  int rank, size, data;
  double adata[SIZE] = {0};
  MPI_Status status;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (rank == 0) {
    //Init Array
    for (int i = 0; i < SIZE; i++) {
      adata[i] = i;
    }
    printf("I have rank %d and I initialized Array :\n", rank);
    print_array(adata);
    MPI_Send(adata, SIZE, MPI_INT, 1,95, MPI_COMM_WORLD);
    MPI_Recv(adata, SIZE, MPI_INT, 1,95, MPI_COMM_WORLD, &status);
    printf("My new array\n", );
    print_array(adata);
  } else {
    MPI_Recv(adata, SIZE, MPI_INT, 0,95, MPI_COMM_WORLD, &status);
    printf("I have rank %d and I receive array:\n", rank);
    print_array(adata);
  }

  MPI_Finalize();
}
