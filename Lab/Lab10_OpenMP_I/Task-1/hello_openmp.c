#include <stdio.h>
#include <omp.h>

//echo $OMP_NUM_THREADS;
//export OMP_NUM_THREADS=3;
//echo $OMP_NUM_THREADS;
int main(int argc, char** argv) {

//printf("OMP_NUM_THREADS: %s\n", getenv("OMP_NUM_THREADS"));

#pragma omp parallel num_threads(10)
  {
    printf("Bonjour!\n");
  }

  return 0;
}
