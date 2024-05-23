#include <stdio.h>
#include <omp.h>

typedef struct DataForThread {
  double A;
  double B;
  int x;
  int y;
} DataForThread_t;

void thread_func(DataForThread_t* p) {
  printf("This is inside thread_func()!\n");
}

int main(int argc, char** argv) {

  const int nThreads = 2;
  DataForThread_t arr[nThreads];
  /* Prepare data for thread 0. */
  arr[0].A = 3.7;
  arr[0].B = 1.2;
  arr[0].x = 88;
  arr[0].y = 77;
  printf("For nThreead=1\narr[0].A:%lf narr[0].B:%lf arr[0].x:%d arr[0].y:%d", arr[0].A, arr[0].B, arr[0].x, arr[0].y);
  /* Prepare data for thread 1. */
  arr[1].A = 2.2;
  arr[1].B = 8.8;
  arr[1].x = 444;
  arr[1].y = 555;
 //printf("For nThreead=1\narr[1].A:%lf", arr[1].A, ":narr[1].B %lf", arr[1].B, ":%d", arr[1].x , ":%d", arr[1].y);
  printf("For nThreead=1\narr[1].A:%lf narr[1].B:%lf arr[1].x:%d arr[1].y:%d", arr[1].A, arr[1].B, arr[1].x, arr[1].y);


#pragma omp parallel num_threads(nThreads)
  {
    /* OK, now we are inside the parallel block, so this is executed
       by several threads. Get ID of current thread. */
    int id = omp_get_thread_num();
    /* Call thread_func and give it a pointer to arr[id] as input. */
    thread_func(&arr[id]);
  }

  return 0;
}
