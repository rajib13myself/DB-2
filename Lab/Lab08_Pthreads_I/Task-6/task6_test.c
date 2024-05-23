#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <pthread.h>
#include <stdlib.h>
#include <unistd.h> // for getopt

#define MAX_THREADS 2

// Structure for passing arguments to the thread function
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

// Function to check if a number is prime
bool is_prime(int num) {
    if (num <= 1) {
        return false;
    }

    for (int i = 2; i * i <= num; i++) {
        if (num % i == 0) {
            return false;
        }
    }
    return true;
}

// Thread function to count prime numbers in a subrange
void* count_primes_thread(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    data->count = 0;
    for (int i = data->start; i <= data->end; i++) {
        if (is_prime(i)) {
            data->count++;
        }
    }
    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    int M = 1000; // Default value for M
    int N = 2;    // Default value for N

    // Parse command line arguments for M and N
    int opt;
    while ((opt = getopt(argc, argv, "M:N:")) != -1) {
        switch (opt) {
            case 'M':
                M = atoi(optarg);
                break;
            case 'N':
                N = atoi(optarg);
                break;
            default:
                fprintf(stderr, "Usage: %s [-M max_number] [-N num_threads]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    pthread_t threads[MAX_THREADS];
    ThreadData thread_data[MAX_THREADS];

    // Determine the number of threads to create
    int num_threads = (N > MAX_THREADS) ? N : MAX_THREADS;

    // Divide the range [1, M] into equal parts for each thread
    int subrange_size = M / num_threads;
    int remainder = M % num_threads; // Remainder for distributing the work

    // Create threads and assign work
    int current_start = 2; // Start from 2 to skip checking 1
    for (int i = 0; i < num_threads; i++) {
        thread_data[i].start = current_start;
        thread_data[i].end = current_start + subrange_size - 1 + (i < remainder ? 1 : 0);
        current_start = thread_data[i].end + 1;

        pthread_create(&threads[i], NULL, count_primes_thread, (void*)&thread_data[i]);
    }

    // Join threads and accumulate counts
    int total_count = 0;

    double start_time = get_wall_seconds();
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
        total_count += thread_data[i].count;
    }
    double end_time = get_wall_seconds();

    printf("Number of prime numbers from 1 to %d: %d with computation time: %0.6f\n", M, total_count, end_time - start_time);

    return 0;
}

