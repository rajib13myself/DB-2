#include <stdio.h>
#include <pthread.h>
#include <time.h>

const long int N1 = 700000000;
const long int N2 = 100000000;

double get_wall_seconds() {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + ts.tv_nsec / 1e9;
}

void* the_thread_func(void* arg) {
  long int i;
  long int sum = 0;
  
  double start_time = get_wall_seconds(); // Start timing the thread function

  for(i = 0; i < N2; i++)
    sum += 7;

  double end_time = get_wall_seconds(); // End timing the thread function
  
  printf("Thread function execution time: %.6f seconds\n", end_time - start_time);

  /* OK, now we have computed sum. Now copy the result to the location given by arg. */
  printf("Value of sum at the_thread_func before assign: %ld\n", sum);
  long int * resultPtr;
  resultPtr = (long int *)arg;
  *resultPtr = sum;
  return NULL;
}

int main() {
  printf("This is the main() function starting.\n");

  long int thread_result_value = 0;

  /* Start thread. */
  pthread_t thread;
  printf("the main() function now calling pthread_create().\n");
  pthread_create(&thread, NULL, the_thread_func, &thread_result_value);

  printf("This is the main() function after pthread_create()\n");

  long int i;
  long int sum = 0;
  double start_time = get_wall_seconds(); // Start timing the main function
  
  for(i = 0; i < N1; i++)
    sum += 7;

  double end_time = get_wall_seconds(); //End timing the main function

  printf("The main() function execution time: %0.6f seconds\n", end_time - start_time);

  /* Wait for thread to finish. */
  printf("the main() function now calling pthread_join().\n");
  pthread_join(thread, NULL);

  printf("sum computed by main() : %ld\n", sum);
  printf("sum computed by thread : %ld\n", thread_result_value);
  long int totalSum = sum + thread_result_value;
  printf("totalSum : %ld\n", totalSum);

  return 0;
}
