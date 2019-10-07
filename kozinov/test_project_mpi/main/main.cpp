#include "mpi.h"
#include <stdio.h>
#include "add.h"


int main(int argc, char ** argv)
{
  MPI_Init(&argc, &argv);
  printf("hello!\n");
  printf("%d\n", add(1, 2));
  MPI_Finalize();
  return 0;
}