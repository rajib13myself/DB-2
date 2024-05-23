#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <pthread.h>
#include <stdlib.h>

#define MAX_THREADS 2

//Structure for passing arguments to the thread function
typedef struct {
  int start;
  int end;
  int count;
} ThreadData;


double get_wall_seconds() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec / 1e9;
}

//Fucntion check the number is prime.
bool is_prime(int num) {
  if (num <= 1) {
    return false;
  }

  for (int i = 2; i * i <=num; i++) {
     if( num % i == 0) {
        return false;
     }
  }
  return true;
}

//Thread Function count the prime numbers in the range of 1 to M
void* count_primes_thread(void* arg) {
  ThreadData* data = (ThreadData*)arg;

  data->count = 0; //Initialize the count variable
  for (int i = data->start; i <= data->end; i++) {
     if(is_prime(i)) {
	data->count++;
     }
  }
  pthread_exit(NULL);
}

int main() {
  int M;
  printf("Please Enter a limit value: ");
  scanf("%d", &M);
  
  pthread_t threads[MAX_THREADS];
  ThreadData thread_data[MAX_THREADS];
  
  //Divide the range[1,M] into two equal parts
  int mid = M / 2;
  
  //Create Thread and Assign work
  for(int i = 0; i < MAX_THREADS; i++) {
     if(i == 0) {
	thread_data[i].start = 2;
	thread_data[i].end = mid;
     } else {
	thread_data[i].start = mid + 1;
	thread_data[i].end = M;
     }
     pthread_create(&threads[i], NULL, count_primes_thread, (void*)&thread_data[i]);
  }
  //Join Threads and accumulate counts
  int total_count = 0;

  double start_time = get_wall_seconds();
  for(int i = 0; i < MAX_THREADS; i++) {
	pthread_join(threads[i], NULL);
	total_count += thread_data[i].count;
  }
  double end_time = get_wall_seconds();

  printf("Number of prime numbers from 1 to %d: %d with computation time: %0.6f\n", M, total_count, end_time-start_time);
  
  return 0;

}
