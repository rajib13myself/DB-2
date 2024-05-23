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
  int thread_index;
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

//Thread Function count the prime numbers in the range of 1 to M from Task-4
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

//Thread Function for number of threads from Task-5
void* thread_function(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    printf("Thread %d\n", data->thread_index);
    pthread_exit(NULL);
}

int main() {
  int M;
  int N;
  printf("Please Enter a limit value (M) and number of Threads (N) : ");
  scanf("%d %d", &M, &N);

  if (M > 0) {
    // Case M
    pthread_t* threads = malloc(MAX_THREADS * sizeof(pthread_t));
    ThreadData* thread_data = malloc(MAX_THREADS * sizeof(ThreadData));
    if (threads == NULL || thread_data == NULL) {
      printf("Memory allocation failed.\n");
      return 1;
    }

    // Divide the range[1,M] into two equal parts
    int mid = M / 2;

    // Create Thread and Assign work
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

    // Join Threads and accumulate counts
    int total_count = 0;
    double start_time = get_wall_seconds();
    for(int i = 0; i < MAX_THREADS; i++) {
      pthread_join(threads[i], NULL);
      total_count += thread_data[i].count;
    }
    double end_time = get_wall_seconds();

    printf("Number of prime numbers from 1 to %d: %d with computation time: %0.6f\n", M, total_count, end_time-start_time);

    free(threads);
    free(thread_data);
  } else if (N > 0) {
    // Case N
    pthread_t* threads = malloc(N * sizeof(pthread_t));
    ThreadData* thread_data = malloc(N * sizeof(ThreadData));
    if (threads == NULL || thread_data == NULL) {
      printf("Memory allocation failed.\n");
      return 1;
    }

    // Create threads
    for (int i = 0; i < N; i++) {
      thread_data[i].thread_index = i;
      pthread_create(&threads[i], NULL, thread_function, (void*)&thread_data[i]);
    }

    // Join threads
    for (int i = 0; i < N; i++) {
      pthread_join(threads[i], NULL);
    }

    free(threads);
    free(thread_data);
  } else {
    printf("Invalid input. Both M and N must be positive integers.\n");
    return 1;
  }

  return 0;
}




