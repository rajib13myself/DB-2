#include <stdio.h>
#include <omp.h>

void thread_func() {
  //int n;
// Get the thread ID number
    int tid = omp_get_thread_num();
    
    // Get the total number of active threads
    int num_threads = omp_get_num_threads();
    
    printf("Thread ID: %d, Total threads: %d\n", tid, num_threads);
   printf("This is inside thread_func() for thread %d!\n", tid);
 // printf("This is inside thread_func(n)!\n");
}

int main(int argc, char** argv) {


#pragma omp parallel num_threads(4)
  {
    thread_func();
  }

  return 0;
}
